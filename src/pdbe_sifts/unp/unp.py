"""
Filters information from UniProt XML and stores it in a dictionary.
This allows for faster access to a small subset of required information
from a large XML file.
"""

import os
import pickle
import random
from collections.abc import Mapping
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any

import funcy
import requests
from funcy.calc import memoize
from lxml import etree as ElementTree
from lxml.etree import XMLSyntaxError

from pdbe_sifts.base.exceptions import ObsoleteUniProtError
from pdbe_sifts.base.log import logger
from pdbe_sifts.base.paths import uniprot_cache_dir as get_uniprot_cache_dir
from pdbe_sifts.base.utils import fetch_uniprot_file

COLORS = {
    "white": 0,
    "red": 31,
    "green": 32,
    "yellow": 33,
    "blue": 34,
}


@memoize
def get_unp_object(acc: str) -> "UNP | None":
    """Return a memoised :class:`UNP` instance for the given accession.

    Handles obsolete accessions gracefully by returning ``None`` and
    skipping the cache entry.  Other exceptions are propagated after also
    skipping the cache so that callers can retry.

    Args:
        acc (str): A UniProt accession (e.g. ``"P12345"`` or
            ``"P12345-2"`` for an isoform).

    Returns:
        UNP | None: The populated :class:`UNP` object, or ``None`` if the
        accession is obsolete.

    Raises:
        Exception: Any exception raised during :class:`UNP` construction
            other than :exc:`~pdbe_sifts.base.exceptions.ObsoleteUniProtError`.
    """
    unp = None
    try:
        unp = UNP(acc)
        logger.info(f"Got UniProt for {acc}")
    except ObsoleteUniProtError:
        logger.warning(f"Obsolete UniProt accession: {acc}. Will be ignored")
    except Exception:
        logger.error(f"Could not get UniProt for {acc}", exc_info=True)
        memoize.skip()
        raise

    if not unp or os.getenv("SIFTS_NO_CACHE_ALL"):
        memoize.skip()
    return unp


def get_annotation_score(accession: str) -> float | None:
    """
    Fetch UniProt annotation score one accession.
    """
    if not accession:
        return []
    url = f"https://rest.uniprot.org/uniprotkb/{accession}?fields=annotation_score"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    data = r.json()
    annotation_score = data.get("annotationScore")
    return annotation_score


def fetch_uniprot_batch(accessions: list[str]) -> list[dict]:
    """
    Fetch UniProt info for up to 100 accessions.
    Returns list of dicts.
    """
    if not accessions:
        return []

    query = "+OR+".join(f"accession:{acc}" for acc in accessions)
    url = (
        "https://rest.uniprot.org/uniprotkb/stream"
        f"?query=({query})"
        "&fields=annotation_score,xref_pdb"
    )

    r = requests.get(url, timeout=30)
    r.raise_for_status()
    data = r.json()

    rows = []

    for e in data.get("results", []):
        accession = e.get("primaryAccession")

        entry_type = e.get("entryType", "")
        provenance = "Swiss-Prot" if "Swiss-Prot" in entry_type else "TrEMBL"

        xrefs = e.get("uniProtKBCrossReferences", [])
        pdb_xref = sum(1 for x in xrefs if x.get("database") == "PDB")

        annotation_score = e.get("annotationScore")

        rows.append(
            {
                "accession": accession,
                "provenance": provenance,
                "pdb_xref": pdb_xref,
                "annotation_score": annotation_score,
            }
        )

    return rows


def fetch_accessions(
    all_accessions: list[str], max_workers: int = 4, batch_size: int = 100
) -> list[dict]:
    """Fetch UniProt metadata for a list of accessions in parallel batches.

    Splits ``all_accessions`` into chunks of at most ``batch_size`` and
    submits each chunk to :func:`fetch_uniprot_batch` via a thread pool.

    Args:
        all_accessions (list[str]): UniProt accessions to look up.
        max_workers (int): Maximum number of concurrent threads.
            Defaults to 4.
        batch_size (int): Number of accessions per API request.
            Defaults to 100.

    Returns:
        list[dict]: Combined results from all batches.  Each dict contains
        ``accession``, ``provenance``, ``pdb_xref``, and
        ``annotation_score`` keys.
    """
    batches = [
        all_accessions[i : i + batch_size]
        for i in range(0, len(all_accessions), batch_size)
    ]

    results = []

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(fetch_uniprot_batch, b) for b in batches]
        for f in as_completed(futures):
            results.extend(f.result())

    return results


def colored(string: str, color: str = "green", bold: bool = False) -> str:
    """Wrap a string in ANSI color escape codes.

    Args:
        string (str): The text to colorise.
        color (str): One of ``"white"``, ``"red"``, ``"green"``,
            ``"yellow"``, or ``"blue"``.  Defaults to ``"green"``.
        bold (bool): When ``True``, prepend the bold escape code.
            Defaults to ``False``.

    Returns:
        str: The input string wrapped in the appropriate ANSI codes.
    """
    if bold:
        b = "\033[1m"
    else:
        b = ""

    tone = f"\033[{COLORS[color]}m"
    reset = "\033[0m"

    return b + tone + str(string) + reset


def in_blue(text: str) -> str:
    """Return ``text`` wrapped in the ANSI blue color code.

    Args:
        text (str): The string to colorise.

    Returns:
        str: The colorised string.
    """
    return colored(text, "blue")


def in_red(text: str) -> str:
    """Return ``text`` wrapped in the ANSI red color code.

    Args:
        text (str): The string to colorise.

    Returns:
        str: The colorised string.
    """
    return colored(text, "red")


def in_green(text: str) -> str:
    """Return ``text`` wrapped in the ANSI green color code.

    Args:
        text (str): The string to colorise.

    Returns:
        str: The colorised string.
    """
    return colored(text, "green")


def in_yellow(text: str) -> str:
    """Return ``text`` wrapped in the ANSI yellow color code.

    Args:
        text (str): The string to colorise.

    Returns:
        str: The colorised string.
    """
    return colored(text, "yellow")


def in_white(text: str) -> str:
    """Return ``text`` wrapped in the ANSI white color code.

    Args:
        text (str): The string to colorise.

    Returns:
        str: The colorised string.
    """
    return colored(text, "white")


def get_item(item: Any) -> str | None:
    """Return the first element of ``item`` as a string, or ``None`` if empty.

    Args:
        item: Any indexable sequence.

    Returns:
        str | None: ``str(item[0])`` if the sequence is non-empty, else
        ``None``.
    """
    try:
        return str(item[0])
    except IndexError:
        return


class UNP:
    """Represents a single UniProt entry with parsed sequence and annotation data.

    Data is sourced from the UniProt XML API.  On first access the XML is
    fetched and cached as a pickle file; subsequent accesses load from that
    cache unless the ``SIFTS_NO_CACHE_ALL`` environment variable is set.

    Attributes:
        accession (str): Primary UniProt accession.
        ad_dbref_auto_acc (str): The accession as originally supplied to
            :meth:`__init__` (may include an isoform suffix such as ``-2``).
        sequence (str): Canonical one-letter-code protein sequence.
        keywords (list[str] | None): Upper-cased UniProt keyword list.
        signal (tuple[int | None, int | None] | None): Signal-peptide
            ``(begin, end)`` positions, or ``None`` when absent.
        chain (tuple[int, int]): Chain ``(begin, end)`` positions.
        molName (str): Recommended protein full name (upper-cased).
        synonyms (list[str]): Short/alternative names from the XML.
        longName (str): UniProt entry name (e.g. ``"SPIKE_SARS2"``).
        organism (list[str]): Scientific and common organism names.
        taxonomy (list[str]): NCBI Taxonomy IDs.
        date_created (str): Entry creation date string.
        date_annotated (list[str]): ``[modified_date, version]``.
        date_seq_update (list[str]): ``[modified_date, version]`` for the
            sequence.
        seq_length (int): Sequence length in residues.
        seq_molWeight (str): Molecular weight from the XML.
        seq_checksum (str): CRC64 checksum of the sequence.
        comments (Mapping[str, list[str]]): Parsed ``<comment>`` blocks
            keyed by comment type.
        dbentry_id (str | None): Internal database entry ID, if set.
        dbreferences (Mapping[str, list[str]]): Cross-references keyed by
            database type.
        features (Mapping[str, list[dict]]): Feature annotations keyed by
            feature type (e.g. ``"sequence variant"``, ``"modified residue"``).
        dataset (str): ``"Swiss-Prot"`` or ``"TrEMBL"``.
        evidences (list[dict]): Parsed evidence blocks.
        secondary_accessions (list[str]): Secondary accessions from the XML.
        seq_isoforms (Mapping[str, str]): Canonical and alternative isoform
            sequences keyed by accession string.
        isoforms (Mapping[str, Any]): Isoform metadata from
            ``<comment type="alternative products">``.
    """

    accession: str
    ad_dbref_auto_acc: str
    sequence: str
    keywords: str
    signal: str
    chain: str
    molName: str
    synonyms: list[str] = []
    longName: str
    organism: str
    taxonomy: str
    date_created: str
    date_annotated: str
    date_seq_update: str
    seq_length: str
    seq_molWeight: str
    seq_checksum: str
    comments: Mapping[str, str] = {}
    dbentry_id: str = None
    dbreferences: Mapping[str, list[str]] = {}
    features: Mapping[str, str] = {}
    dataset: str
    evidences: list[str] = []
    secondary_accessions: list[str] = []

    def _load_uniprot_data(self, accession: str) -> None:
        """Load the data for the accession from the uniprot XML file

        Will try to load from a pickle file if it exists. If not, will fetch from uniprot and cache it.
        if env var SIFTS_NO_CACHE_ALL is set to True, will not cache the data nor load from cache.

        Args:
            accession (str): The accession to load

        """
        cache = bool(os.getenv("SIFTS_NO_CACHE_ALL", True))

        pickle_file_path = os.path.join(
            get_uniprot_cache_dir(accession), f"{accession}.pkl"
        )

        loaded = False
        if cache and Path(pickle_file_path).exists():
            loaded = self._load_from_pickle(accession, pickle_file_path)
            if loaded:
                return

        self._load_from_uniprot_api(accession)

        if cache:
            Path(pickle_file_path).parent.mkdir(parents=True, exist_ok=True)
            with open(pickle_file_path, "wb") as f:
                pickle.dump(self.__dict__, f, pickle.HIGHEST_PROTOCOL)

    def _load_from_pickle(self, accession: str, pickle_file_path: str) -> bool:
        """Attempt to populate this instance from a cached pickle file.

        If the pickle cannot be loaded (e.g. it is corrupt or from an
        incompatible version), the file is deleted so that a fresh fetch
        is attempted on the next call.

        Args:
            accession (str): The UniProt accession, used only for logging.
            pickle_file_path (str): Absolute path to the pickle cache file.

        Returns:
            bool: ``True`` if the pickle was loaded successfully, ``False``
            otherwise.
        """
        with open(pickle_file_path, "rb") as f:
            try:
                self.__dict__ = pickle.load(f)
                logger.info(f"Loaded UNP pickle from {pickle_file_path}")
                return True

            except Exception:
                logger.warning(
                    f"Could not load pickle file for {accession}. Will delete it and fetch from uniprot"
                )
                os.unlink(pickle_file_path)
        return False

    @funcy.retry(3, errors=XMLSyntaxError, timeout=0.2)
    def parse_xml(self, xml_file: str) -> Any:
        """Parse the XML file and populate the fields

        Retries 3 times if there is an XMLSyntaxError to avoid
        race condition flushing the file to disk

        Args:
            xml_file (str): Path to the XML file
        """
        return ElementTree.parse(xml_file).getroot()

    def _load_from_uniprot_api(self, accession: str) -> dict:
        """Fetch the UniProt XML for ``accession`` and populate this instance.

        Downloads the XML via :func:`~pdbe_sifts.base.utils.fetch_uniprot_file`,
        parses it with :meth:`parse_xml`, and delegates field extraction to
        :meth:`_populate_fields`.

        Args:
            accession (str): The UniProt accession to fetch.

        Raises:
            ValueError: If the downloaded XML contains no accession elements.
        """
        xml_file = fetch_uniprot_file(accession, "xml")
        doc = self.parse_xml(xml_file)

        accessions = list(
            map(str, doc[0].xpath(".//*[name()='accession']/text()"))
        )
        if not accessions:
            raise ValueError(
                f"Invalid file for {accession}. No accessions found."
            )

        self.accession = accession
        accessions = list(
            map(str, doc[0].xpath(".//*[name()='accession']/text()"))
        )
        if not accessions:
            raise ValueError(
                f"Invalid file for {accession}. No accessions found."
            )

        self.accession = accession
        self._populate_fields(doc)

    def __init__(self, accession: str) -> None:
        """Initialise and fully populate a :class:`UNP` instance.

        Strips any isoform suffix (``-N``) before fetching data so that the
        canonical entry is always retrieved, while preserving the original
        accession string in :attr:`ad_dbref_auto_acc`.

        Args:
            accession (str): A UniProt accession, optionally with an isoform
                suffix (e.g. ``"P12345-2"``).
        """
        self.seq_isoforms: Mapping[str, Any] = {}
        self.isoforms: Mapping[str, Any] = {}
        self.features: Mapping[str, str] = {}
        self.accession = None
        self.ad_dbref_auto_acc = accession
        if "-" in accession:
            accession = accession.split("-")[0]
        self._load_uniprot_data(accession)

    def _populate_fields(self, doc: Any) -> None:
        """Extract and store all fields from a parsed UniProt XML document.

        Populates accession, sequence, keywords, signal peptide, chain,
        protein names, organism, taxonomy, entry dates, comments,
        cross-references, features, and evidences from ``doc``.

        Args:
            doc: The root element returned by :meth:`parse_xml`.
        """
        uniprot = doc[0]
        self.secondary_accessions = list(
            map(str, uniprot.xpath(".//*[name()='accession']/text()"))
        )
        self.accession = str(self.secondary_accessions[0])
        del self.secondary_accessions[0]

        about_sequence = uniprot.xpath(
            ".//*[name()='sequence' and ./@length!='']"
        )[0]
        self.sequence = about_sequence.xpath("./text()")[0].replace("\n", "")

        self.date_seq_update = list(
            map(
                str,
                (
                    about_sequence.xpath("./@modified")[0],
                    about_sequence.xpath("./@version")[0],
                ),
            )
        )
        self.seq_length = int(about_sequence.xpath("./@length")[0])
        self.seq_molWeight = str(about_sequence.xpath("./@mass")[0])
        self.seq_checksum = str(about_sequence.xpath("./@checksum")[0])

        self.keywords = [
            kw.upper() for kw in uniprot.xpath(".//*[name()='keyword']/text()")
        ]

        if len(self.keywords) < 1:
            self.keywords = None

        self.signal = self._extract_signal(
            uniprot.xpath(".//*[name()='feature'][@type='signal peptide']")
        )

        chainIndices = uniprot.xpath(".//*[name()='feature'][@type='chain']")

        if chainIndices != [] and (
            chainIndices[0].xpath(".//*[name()='begin']/@status") != []
            and chainIndices[0].xpath(".//*[name()='end']/@status") != []
        ):
            self.chain = (
                int(chainIndices[0].xpath(".//*[name()='begin']/@position")[0]),
                int(chainIndices[0].xpath(".//*[name()='end']/@position")[0]),
            )

        self.molName = [
            name.upper()
            for name in uniprot.xpath(".//*[name()='protein']//text()")
        ]

        recommendedNames = uniprot.xpath(".//*[name()='recommendedName']")
        alternativeNames = uniprot.xpath(".//*[name()='alternativeName']")
        if recommendedNames:
            self.synonyms = recommendedNames[0].xpath(
                ".//*[name()!='ecNumber' and name()!='fullName']/text()"
            )
        for alt in alternativeNames:
            self.synonyms += alt.xpath(".//*[name()!='ecNumber']/text()")

        self.molName = self.molName[0]
        self.longName = str(uniprot.xpath(".//*[name()='name']/text()")[0])
        self._extract_taxonomy(uniprot.xpath(".//*[name()='organism']"))

        about_entry = doc.xpath(".//*[name()='entry']")[0]
        self.dataset = str(about_entry.xpath("./@dataset")[0])
        self.date_created = str(about_entry.xpath("./@created")[0])
        self.date_annotated = list(
            map(
                str,
                (
                    about_entry.xpath("./@modified")[0],
                    about_entry.xpath("./@version")[0],
                ),
            )
        )

        about_comments = uniprot.xpath(".//*[name()='comment']")
        self._extract_comments(about_comments)

        self.dbreferences = {}
        about_dbreference = uniprot.xpath(".//*[name()='dbReference']")
        for ref in about_dbreference:
            ref_type = ref.xpath("./@type")[0]
            self.dbreferences.setdefault(ref_type, []).append(
                str(ref.xpath("./@id")[0])
            )

        about_features = uniprot.xpath(".//*[name()='feature']")
        self._extract_features(about_features)

        about_evidences = uniprot.xpath(".//*[name()='evidence']")
        for ev in about_evidences:
            evidence = {}
            for k in ["key", "category", "type", "attribute", "date"]:
                evidence[k] = get_item(ev.xpath(f"./@{k}"))
            self.evidences.append(evidence)

        self.seq_isoforms[self.accession] = self.sequence

    def _extract_taxonomy(self, about_organism: list) -> None:
        """Parse organism and taxonomy data from the ``<organism>`` XML elements.

        Sets :attr:`organism` to a list of scientific and common names, and
        :attr:`taxonomy` to a list of NCBI Taxonomy ID strings.

        Args:
            about_organism (list): XPath result for ``.//*[name()='organism']``
                elements.
        """
        self.organism = about_organism[0].xpath(
            ".//*[name()='name'][@type='scientific']/text()"
        ) + about_organism[0].xpath(
            ".//*[name()='name'][@type='common']/text()"
        )

        self.organism = list(map(str, self.organism))
        self.taxonomy = [
            str(tax.xpath("./@id")[0])
            for tax in about_organism[0].xpath(
                ".//*[name()='dbReference'][@type='NCBI Taxonomy']"
            )
        ]

    def _extract_features(self, about_features: list) -> None:
        """Parse feature annotations from ``<feature>`` XML elements.

        Handles ``"sequence variant"``, ``"modified residue"``,
        ``"mutagenesis site"``, and ``"splice variant"`` feature types,
        storing each as a dict in the corresponding list under
        :attr:`features`.

        Args:
            about_features (list): XPath result for
                ``.//*[name()='feature']`` elements.
        """
        for feat in about_features:
            feat_type = feat.xpath("./@type")[0]

            if feat_type not in list(self.features.keys()):
                self.features[feat_type] = []

            if feat_type in [
                "sequence variant",
                "modified residue",
                "mutagenesis site",
                "splice variant",
            ]:
                if feat.xpath("./@id"):
                    id_value = str(feat.xpath("./@id")[0])
                else:  # check out P26039
                    id_value = str(random.randint(0, 1000))

                if feat.xpath(".//*[name()='original']/text()") == []:
                    original = ""
                    variation = ""
                else:
                    original = str(
                        feat.xpath(".//*[name()='original']/text()")[0]
                    )
                    variation = str(
                        feat.xpath(".//*[name()='variation']/text()")[0]
                    )

                if feat.xpath("./@evidence") == []:
                    evidence = None
                else:
                    evidence = str(feat.xpath("./@evidence")[0])

                location = feat.xpath(".//*[name()='location']")[0]

                if location.xpath(".//*[name()='position']") == []:
                    begin = location.xpath(".//*[name()='begin']/@position")[0]
                    end = location.xpath(".//*[name()='end']/@position")[0]
                else:
                    begin = location.xpath(".//*[name()='position']/@position")[
                        0
                    ]
                    end = begin

                if feat_type in [
                    "sequence variant",
                    "mutagenesis site",
                    "splice variant",
                ]:
                    self.features[feat_type].append(
                        {
                            "id": id_value,
                            "evidence": evidence,
                            "original": original,
                            "variation": variation,
                            "begin": int(begin),
                            "end": int(end),
                        }
                    )
                elif feat_type == "modified residue":
                    self.features[feat_type].append(
                        {
                            "id": id_value,
                            "evidence": evidence,
                            "position": int(begin),
                        }
                    )

    def _extract_comments(self, about_comments: list) -> None:
        """Parse comment blocks from ``<comment>`` XML elements.

        Populates :attr:`comments` and, for ``"alternative products"``
        comments, also populates :attr:`isoforms` with splice-variant
        sequence references.

        Args:
            about_comments (list): XPath result for
                ``.//*[name()='comment']`` elements.
        """
        for cc in about_comments:
            cc_type = cc.xpath("./@type")[0]

            if cc_type not in list(self.comments.keys()):
                self.comments[cc_type] = []

            if cc_type == "alternative products":
                for iso in cc.xpath(".//*[name()='isoform']"):
                    isoid = iso.xpath(".//*[name()='id']")[0]
                    isoid = isoid.xpath("./text()")[0]
                    spliceids = iso.xpath(".//*[name()='sequence']/@ref")

                    if spliceids:
                        spliceids = spliceids[0]
                        spliceids = spliceids.split(" ")

                    self.isoforms[str(isoid)] = spliceids

            if cc.xpath(".//*[name()='text']/text()") != []:
                self.comments[cc_type].append(
                    str(cc.xpath(".//*[name()='text']/text()")[0])
                )
            if cc.xpath(".//*[name()='link']/@uri") != []:
                self.comments[cc_type].append(
                    str(cc.xpath(".//*[name()='link']/@uri")[0])
                )
            if cc.xpath(".//*[name()='subcellularLocation']") != []:
                self.comments[cc_type].append(
                    str(
                        cc.xpath(".//*[name()='subcellularLocation']")[0].xpath(
                            ".//*[name()='location']/text()"
                        )[0]
                    )
                )

    def _extract_signal(
        self, signal_peptide: list
    ) -> tuple[int | None, int | None] | None:
        """Parse the signal-peptide location from the relevant ``<feature>`` element.

        Handles the cases where begin or end positions are marked as
        ``"unknown"`` in the XML, and where a single ``<position>`` element
        is used instead of separate ``<begin>``/``<end>`` elements.

        Args:
            signal_peptide (list): XPath result for
                ``.//*[name()='feature'][@type='signal peptide']``.

        Returns:
            tuple[int | None, int | None] | None: ``(begin, end)`` residue
            positions (1-based), where either value may be ``None`` if the
            position is unknown.  Returns ``None`` if no signal peptide
            feature is present.
        """
        signal = None
        if signal_peptide:
            signal_peptide = signal_peptide[0]

            if (
                signal_peptide.xpath(".//*[name()='begin']/@status")
                and signal_peptide.xpath(".//*[name()='begin']/@status")[0]
                == "unknown"
            ):
                begin = None
            # if begin == end then they sometimes use position/position
            elif signal_peptide.xpath(".//*[name()='position']/@position"):
                begin = int(
                    signal_peptide.xpath(".//*[name()='position']/@position")[0]
                )
            else:
                begin = int(
                    signal_peptide.xpath(".//*[name()='begin']/@position")[0]
                )

            if (
                signal_peptide.xpath(".//*[name()='end']/@status")
                and signal_peptide.xpath(".//*[name()='end']/@status")[0]
                == "unknown"
            ):
                end = None
            elif signal_peptide.xpath(".//*[name()='position']/@position"):
                end = begin
            else:
                end = int(
                    signal_peptide.xpath(".//*[name()='end']/@position")[0]
                )

            signal = (begin, end)

        return signal

    def getAllIsoforms(self) -> dict[str, str]:
        """Return all known isoform sequences keyed by isoform accession.

        If isoforms have not yet been computed, applies each known splice
        variant from :attr:`isoforms` to generate the sequences.

        Returns:
            dict[str, str]: Mapping of isoform accession to one-letter-code
            sequence.
        """
        if len(self.seq_isoforms) > 1:
            return self.seq_isoforms

        for iso, variants in list(self.isoforms.items()):
            # Don't follow isoforms with a different accession
            if iso.split("-")[0] == self.accession:
                self.seq_isoforms[iso] = self.applySpliceVariant(variants)

        return self.seq_isoforms

    def applySpliceVariant(self, spliceids: str | list) -> str | None:
        """Apply one or more splice variant IDs to the canonical sequence.

        Variants are applied in reverse order of their start position so
        that earlier edits do not shift the indices of later ones.

        Args:
            spliceids (str | list): A single splice variant ID string or a
                list of IDs to apply simultaneously.

        Returns:
            str | None: The modified sequence, or ``None`` if no splice
            variants are present in :attr:`features`.
        """
        newSeq = None
        splices = []

        if isinstance(spliceids, str):
            spliceids = [spliceids]

        if self.hasSpliceVariants():
            for splice in spliceids:
                for var in self.features["splice variant"]:
                    if var["id"] == splice:
                        splices.append(var)
                        break

        splices = sorted(splices, key=lambda x: x["begin"], reverse=True)

        newSeq = self.sequence

        for splice in splices:
            newSeq = (
                newSeq[: splice["begin"] - 1]
                + splice["variation"]
                + newSeq[splice["end"] :]
            )

        return newSeq

    def getAllSpliceVariants(self) -> list:
        """Return the canonical sequence plus all individual splice-variant sequences.

        Returns:
            list[str]: The canonical sequence followed by the result of
            applying each splice variant individually.  Returns an empty
            list if no splice variants exist.
        """
        splices = []

        if self.hasSpliceVariants():
            splices.append(self.sequence)

            # splices.append(
            #     self.applySpliceVariant(
            #         [x["id"] for x in self.features["splice variant"]]
            #     )
            # )

            for var in self.features["splice variant"]:
                splices.append(self.applySpliceVariant(var["id"]))

        return splices

    def applySequenceVariant(self, variant: str) -> str | None:
        """Apply a single sequence-variant substitution to the canonical sequence.

        Args:
            variant (str): The variant ID (e.g. ``"VAR_012345"``).

        Returns:
            str | None: The modified sequence, or ``None`` if the variant ID
            is not found or no sequence variants exist.
        """
        newSeq = None
        if self.hasVariants():
            for var in self.features["sequence variant"]:
                if var["id"] == variant:
                    newSeq = (
                        self.sequence[: var["begin"] - 1]
                        + var["variation"]
                        + self.sequence[var["end"] :]
                    )
                    break

        return newSeq

    def getAllSequenceVariants(self) -> list:
        """Return the canonical sequence plus all individual sequence-variant sequences.

        Returns:
            list[str | None]: The canonical sequence followed by the result of
            applying each sequence variant individually.  Returns an empty
            list if no sequence variants exist.
        """
        variants = []

        if self.hasVariants():
            variants.append(self.sequence)

            for var in self.features["sequence variant"]:
                variants.append(self.applySequenceVariant(var["id"]))

        return variants

    def isReviewed(self) -> bool:
        """Return whether this entry is a reviewed Swiss-Prot entry.

        Returns:
            bool: ``True`` if :attr:`dataset` is ``"Swiss-Prot"``.
        """
        return self.dataset == "Swiss-Prot"

    def getOrganismScientific(self) -> list[str] | None:
        """Return the organism name(s) or ``None`` when absent.

        Returns:
            list[str] | None: Scientific and common organism name strings,
            or ``None`` if the list is empty.
        """
        return self.organism if self.organism else None

    def getTaxId(self) -> list[str] | None:
        """Return the NCBI Taxonomy IDs or ``None`` when absent.

        Returns:
            list[str] | None: NCBI Taxonomy ID strings, or ``None`` if the
            list is empty.
        """
        return self.taxonomy if self.taxonomy else None

    def getUNPPrimaryName(self) -> str | None:
        """Return the recommended protein full name or ``None`` when absent.

        Returns:
            str | None: The value of :attr:`molName`, or ``None`` if unset.
        """
        return self.molName if self.molName else None

    def hasSecondaryName(self) -> bool:
        """Return whether any alternative/short names are present.

        Returns:
            bool: ``True`` if :attr:`synonyms` is non-empty.
        """
        return bool(self.synonyms)

    def getUNPSecondaryName(self) -> list[str]:
        """Return alternative/short protein names.

        Returns:
            list[str]: The :attr:`synonyms` list, or an empty list if none
            are recorded.
        """
        return self.synonyms if self.synonyms else []

    def hasShortName(self) -> bool:
        """Return whether any synonyms (used as short names) are present.

        Returns:
            bool: ``True`` if :attr:`synonyms` is non-empty.
        """
        return bool(self.synonyms)

    def getUNPShortName(self) -> str | bool:
        """Return the UniProt entry name (used as a short name) or ``False``.

        Returns:
            str | bool: :attr:`longName` if set, otherwise ``False``.
        """
        return self.longName if self.longName else False

    def hasSpliceVariants(self) -> bool:
        """Return whether any splice variants are annotated.

        Returns:
            bool: ``True`` if :attr:`features` contains a
            ``"splice variant"`` key.
        """
        return bool(
            self.features and "splice variant" in list(self.features.keys())
        )

    def hasVariants(self) -> bool:
        """Return whether any sequence variants are annotated.

        Returns:
            bool: ``True`` if :attr:`features` contains a
            ``"sequence variant"`` key.
        """
        return bool(
            self.features and "sequence variant" in list(self.features.keys())
        )

    def getSpliceVariants(self) -> list:
        """Return all splice-variant feature dicts.

        Returns:
            list[dict]: Each dict contains ``id``, ``evidence``,
            ``original``, ``variation``, ``begin``, and ``end`` keys.
            Returns an empty list when no splice variants exist.
        """
        if self.hasSpliceVariants():
            return self.features["splice variant"]

        return []

    def getVariants(self) -> list:
        """Return all sequence-variant feature dicts.

        Returns:
            list[dict]: Each dict contains ``id``, ``evidence``,
            ``original``, ``variation``, ``begin``, and ``end`` keys.
            Returns an empty list when no sequence variants exist.
        """
        if self.hasVariants():
            return self.features["sequence variant"]

        return []

    def getSequenceOneLetter(self) -> str:
        """Return the canonical one-letter-code protein sequence.

        Returns:
            str: The value of :attr:`sequence`.
        """
        return self.sequence

    def getAccession(self) -> str:
        """Return the primary UniProt accession.

        Returns:
            str: The value of :attr:`accession`.
        """
        return self.accession

    def getAccessionLong(self) -> str:
        """Return the UniProt entry name (the ``longName`` field).

        Returns:
            str: The value of :attr:`longName` (e.g. ``"SPIKE_SARS2"``).
        """
        return self.longName

    def hasSecondaryAccessions(self) -> bool:
        """Return whether any secondary accessions are recorded.

        Returns:
            bool: ``True`` if :attr:`secondary_accessions` is non-empty.
        """
        return bool(
            self.secondary_accessions and self.secondary_accessions != []
        )

    def hasModifiedResidues(self) -> bool:
        """Return whether any modified-residue annotations are present.

        Returns:
            bool: ``True`` if :attr:`features` contains a
            ``"modified residue"`` key.
        """
        return bool(
            self.features and "modified residue" in list(self.features.keys())
        )

    def hasMutagenesisSite(self) -> bool:
        """Return whether any mutagenesis-site annotations are present.

        Returns:
            bool: ``True`` if :attr:`features` contains a
            ``"mutagenesis site"`` key.
        """
        return bool(
            self.features and "mutagenesis site" in list(self.features.keys())
        )

    def getFeatures(self, idx: int) -> dict:
        """Return all feature annotations that overlap a given residue position.

        Checks splice variants, sequence variants, modified residues, and
        mutagenesis sites.  For range-based features the position must lie
        within ``[begin, end]``; for point features (modified residues) it
        must equal ``position`` exactly.

        Args:
            idx (int): 1-based residue position.

        Returns:
            dict[str, list[str]]: Mapping of feature type to a list of
            feature IDs overlapping ``idx``.  Keys are only present when at
            least one matching feature is found.
        """
        features = {}

        if self.hasSpliceVariants():
            features["splice variant"] = []

            for var in self.features["splice variant"]:
                if var["begin"] <= idx <= var["end"]:
                    features["splice variant"].append(var["id"])

        if self.hasVariants():
            features["sequence variant"] = []

            for var in self.features["sequence variant"]:
                if var["begin"] <= idx <= var["end"]:
                    features["sequence variant"].append(var["id"])

        if self.hasModifiedResidues():
            features["modified residue"] = []

            for var in self.features["modified residue"]:
                if var["position"] == idx:
                    features["modified residue"].append(var["id"])

        if self.hasMutagenesisSite():
            features["mutagenesis site"] = []

            for var in self.features["mutagenesis site"]:
                if var["begin"] <= idx <= var["end"]:
                    features["mutagenesis site"].append(var["id"])

        return features
