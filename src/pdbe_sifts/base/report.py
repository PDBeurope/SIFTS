import datetime
import os
from enum import Enum

from sqlalchemy import Column, DateTime, Integer, MetaData, String, Table
from sqlalchemy.exc import DatabaseError, IntegrityError
from sqlalchemy.sql.base import Executable

from pdbe_sifts.base.log import logger
from pdbe_sifts.base.utils import create_db_connection, get_exc_context


report_conn: Executable | None = None
report_to: str = "db"
entry_table: Table | None = None
summary_table: Table | None = None


class EntryStatusType(Enum):
    SUCCESS = "success"
    FAILED = "failed"
    SKIPPED = "skipped"
    KNOWN_EXCEPTION = "known_exception"

    def __repr__(self):
        return self.value


def create_summary_table(conn: Executable) -> Table:
    meta = MetaData()

    summary = Table(
        "orc_task_summary",
        meta,
        Column("context", String(20), nullable=False),
        Column("exec_date", DateTime, default=datetime.datetime.utcnow),
        Column("module", String(20)),
        Column("status", String(30), nullable=False),
        Column("processed", Integer),
        Column("failed", Integer),
        Column("known_exceptions", String(100)),
        Column("log_dir", String(100)),
        Column("message", String(2000)),
    )

    summary.create(bind=conn, checkfirst=True)
    return summary


def create_log_table(conn: Executable) -> Table:
    meta = MetaData()

    entry_log = Table(
        "orc_entry_log",
        meta,
        Column("context", String(20), nullable=False),
        Column("exec_date", DateTime, default=datetime.datetime.utcnow),
        Column("entry", String(30)),
        Column("module", String(20)),
        Column("status", String(30), nullable=False),
        Column("exception_type", String(100)),
        Column("exception_message", String(2000)),
    )

    entry_log.create(bind=conn, checkfirst=True)
    return entry_log


def setup_entry_logging():
    global report_conn, report_to, entry_table, summary_table
    db_conn_str = os.getenv("ORC_LOG_DB", None)
    if db_conn_str:
        try:
            report_conn = create_db_connection(db_conn_str)
            entry_table = create_log_table(report_conn)
            summary_table = create_summary_table(report_conn)
        except DatabaseError:
            logger.error(
                "Could not create entry_log table. Logging to stdout", exc_info=True
            )
            report_to = "stdout"
    else:
        logger.info("No db credentials specified for entry logging. Will use stdout")


def log_entry(
    entry: str,
    module: str,
    status: EntryStatusType,
    exception_type=None,
    exception_message=None,
):
    global entry_table
    # Setup if not already done
    if report_to == "db":
        if not report_conn:
            setup_entry_logging()
        if entry_table is None and report_conn:
            entry_table = create_log_table(report_conn)

    if entry_table is None or report_conn is None or report_to == "stdout":
        logger.debug(f"Entry_table: {entry_table}; conn:{report_conn}")
        logger.info(
            f"{module}, {entry}, {status}, {exception_type}, {exception_message}"
        )
        return

    context, exec_date = get_exc_context()

    if status.value == "success":
        return

    values = {
        "context": context,
        "exec_date": datetime.datetime.strptime(exec_date, "%Y-%m-%d %H:%M:%S"),
        "status": status.value,
        "entry": entry,
        "module": module,
        "exception_type": str(exception_type) if exception_type else None,
        "exception_message": str(exception_message) if exception_message else None,
    }
    try:
        stmt = entry_table.insert().values(values)
        report_conn.execute(stmt)

    except IntegrityError:
        stmt = entry_table.update().values(values)
        report_conn.execute(stmt)
    # report_conn.commit()


def send_summary_message(
    module,
    status: EntryStatusType,
    processed: int = None,
    failed: int = None,
    known_exceptions: int = None,
    log_dir: str = None,
    message=None,
):
    global summary_table

    context, exec_date = get_exc_context()
    values = {
        "context": context,
        "exec_date": datetime.datetime.strptime(exec_date, "%Y-%m-%d %H:%M:%S"),
        "module": module,
        "status": status.value,
        "processed": processed,
        "failed": failed,
        "known_exceptions": str(known_exceptions),
        "log_dir": log_dir,
        "message": message,
    }

    if report_to == "db":
        if not report_conn:
            setup_entry_logging()
        if summary_table is None and report_conn:
            summary_table = create_log_table(report_conn)

    if summary_table is None or report_conn is None or report_to == "stdout":
        logger.debug(f"summary_table: {summary_table}; conn:{report_conn}")
        logger.info(values)
        return

    try:
        stmt = summary_table.insert().values(values)
        report_conn.execute(stmt)

    except IntegrityError:
        stmt = summary_table.update().values(values)
        report_conn.execute(stmt)
