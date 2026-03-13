import argparse
import shutil
import subprocess
from pathlib import Path

from pdbe_sifts.base.log import logger
from pdbe_sifts.config import load_config

conf = load_config()


class CopyToExchange:
    def __init__(self, input_dir, output_dir, remote=False) -> None:
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.remote = remote
        if remote:
            self.remote_server = conf.exchange.server
            self.username = conf.exchange.username

    def process_entry(self, entry_id: str):
        sifts_only_file = str(Path(self.input_dir) / f"{entry_id}_sifts_only.cif.gz")
        div_output_dir = f"{self.output_dir}/{entry_id[1:3]}"

        if not Path(sifts_only_file).exists():
            logger.warning(f"File {sifts_only_file} does not exist")
            return

        if self.remote:
            ssh_cmd = (
                f'ssh {self.username}@{self.remote_server} "mkdir -p {div_output_dir}"'
            )
            logger.info(f"Create remote output dir using cmd: {ssh_cmd}")
            subprocess.run(ssh_cmd, check=True, shell=True)

            dest = f"{self.username}@{self.remote_server}:{self.output_dir}"
            cmd = f"rsync -av {sifts_only_file} {dest}"
            logger.info(f"Copy file to remote using cmd: {cmd}")
            subprocess.run(cmd, check=True, shell=True)
        else:
            Path(div_output_dir).mkdir(parents=True, exist_ok=True)
            shutil.copy(sifts_only_file, div_output_dir)

    def copy_changes_list(self):
        if self.remote:
            logger.info("Copying list of changes to remote server")
            dest = f"{self.username}@{self.remote_server}:{self.output_dir}/"
            subprocess.run(
                f"scp {conf.lists.sifts_mapping_changes} {dest}",
                check=True,
                shell=True,
            )
        else:
            logger.info("Copying list of changes to output dir")
            shutil.copy(conf.lists.sifts_mapping_changes, self.output_dir)


def run():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--entry",
        required=True,
        help="PDB entry ID to process",
    )
    parser.add_argument(
        "-i",
        "--input-dir",
        default=conf.location.work.data_entry_dir,
        required=True,
        help="Input base directory",
    )
    parser.add_argument(
        "-o",
        "--output-dir",
        help="Output base directory",
    )
    parser.add_argument(
        "--remote",
        action="store_true",
        default=False,
        help=(
            "Copy directly to exchange server. Assumes output-dir "
            "is on the remote server define in conf.exchange section."
        ),
    )

    args = parser.parse_args()

    output_dir = args.output_dir
    if not output_dir:
        output_dir = conf.sifts_to_mmcif.staging_dir
        if args.remote:
            output_dir = conf.sifts_to_mmcif.exchange_dir

    obj = CopyToExchange(args.input_dir, output_dir, args.remote)
    obj.process_entry(args.entry)


if __name__ == "__main__":
    run()
