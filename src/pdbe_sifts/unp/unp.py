#!/usr/bin/env python3

"""
Filters information from UniProt XML and stores it in a dictionary.
This allows for faster access to a small subset of required information
from a large XML file.
"""

import pickle as pkl
import requests
from pathlib import Path
from typing import Any, Optional

import funcy
from lxml import etree as ElementTree
from lxml.etree import XMLSyntaxError
from concurrent.futures import ThreadPoolExecutor, as_completed

from pdbe_sifts.base.log import logger
from pdbe_sifts.base.utils import fetch_uniprot_file, _uniprot_cache_path

# Fields that are instance-specific (set in __init__) and must NOT be
# overwritten from the pickle.
_PICKLE_SKIP_FIELDS = frozenset({"unp_dir", "ad_dbref_auto_acc"})

COLORS = {
    "white": 0,
    "red": 31,
    "green": 32,
    "yellow": 33,
    "blue": 34,
}

def get_annotation_score(accession):
    """
    Fetch UniProt annotation score one accession.
    """
    if not accession:
        return []
    url = f'https://rest.uniprot.org/uniprotkb/{accession}?fields=annotation_score'
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

        rows.append({
            "accession": accession,
            "provenance": provenance,
            "pdb_xref": pdb_xref,
            "annotation_score": annotation_score,
        })

    return rows


def fetch_accessions(all_accessions, max_workers=4, batch_size=100):
    batches = [
        all_accessions[i:i+batch_size]
        for i in range(0, len(all_accessions), batch_size)
    ]

    results = []

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(fetch_uniprot_batch, b) for b in batches]
        for f in as_completed(futures):
            results.extend(f.result())

    return results


def colored(string, color="green", bold=False):
    if bold:
        b = "\033[1m"
    else:
        b = ""

    tone = "\033[%dm" % COLORS[color]
    reset = "\033[0m"

    return b + tone + str(string) + reset


def in_blue(text):
    return colored(text, "blue")


def in_red(text):
    return colored(text, "red")


def in_green(text):
    return colored(text, "green")


def in_yellow(text):
    return colored(text, "yellow")


def in_white(text):
    return colored(text, "white")


def get_item(item):
    try:
        return str(item[0])
    except IndexError:
        return


class UNP:
    accession: str
    unp_dir: str
    ad_dbref_auto_acc: str
    sequence: str
    keywords: Optional[list[str]]
    signal: Optional[tuple]
    chain: Optional[tuple[int, int]]
    molName: str
    synonyms: list[str]
    longName: str
    organism: list[str]
    taxonomy: list[str]
    date_created: str
    date_annotated: list[str]
    date_seq_update: list[str]
    seq_length: int
    seq_molWeight: str
    seq_checksum: str
    comments: dict[str, list]
    dbentry_id: str
    dbreferences: dict[str, list[str]]
    features: dict[str, list]
    dataset: str
    evidences: list[dict]
    secondary_accessions: list[str]
    annotation_score: Optional[int]

    # ------------------------------------------------------------------
    # Cache helpers
    # ------------------------------------------------------------------

    def _pickle_path(self, accession: str) -> Path:
        """Return the pickle cache path for *accession*."""
        return _uniprot_cache_path(accession, "pickle", self.unp_dir)

    def _load_from_pickle(self, path: Path) -> None:
        """Restore instance state from *path* (pickle file)."""
        with open(path, "rb") as fh:
            data: dict = pkl.load(fh)
        # Restore all cached fields, but never overwrite constructor fields.
        for key, value in data.items():
            if key not in _PICKLE_SKIP_FIELDS:
                setattr(self, key, value)

    def _save_to_pickle(self, path: Path) -> None:
        """Persist current instance state to *path* (pickle file)."""
        path.parent.mkdir(parents=True, exist_ok=True)
        data = {k: v for k, v in self.__dict__.items() if k not in _PICKLE_SKIP_FIELDS}
        tmp = path.with_suffix(".tmp")
        try:
            with open(tmp, "wb") as fh:
                pkl.dump(data, fh, protocol=pkl.HIGHEST_PROTOCOL)
            tmp.replace(path)  # atomic rename
        except Exception:
            logger.warning(f"Could not save pickle for {path.stem}")
            tmp.unlink(missing_ok=True)

    # ------------------------------------------------------------------
    # Data loading
    # ------------------------------------------------------------------

    def _load_uniprot_data(self, accession: str, fetch_annotation_score: bool = False) -> None:
        """Load UNP data, preferring the local cache.

        Load order:
        1. Pickle (``{unp_dir}/{first_letter}/{accession}.pickle``) — fastest,
           zero network traffic.
        2. XML on disk (``{unp_dir}/{first_letter}/{accession}.xml``) — parse
           locally, then save a pickle for future runs.
        3. Download XML from UniProt REST API — save XML + pickle.

        Args:
            accession: Canonical UniProt accession (no isoform suffix).
            fetch_annotation_score: When *True* and the pickle does not already
                contain ``annotation_score``, fetch it from the API and include
                it in the newly written pickle.
        """
        pickle_path = self._pickle_path(accession)

        if pickle_path.exists():
            logger.debug(f"{accession}: loading from pickle cache")
            self._load_from_pickle(pickle_path)
            # annotation_score may already be populated from the pickle.
            return

        # Parse from XML (downloads if not yet cached).
        self._load_from_uniprot_xml(accession)

        # Optionally fetch annotation score while we are about to save the pickle
        # so that future loads need zero network round-trips.
        if fetch_annotation_score and self.annotation_score is None:
            self.annotation_score = get_annotation_score(accession)

        self._save_to_pickle(pickle_path)

    @funcy.retry(3, errors=XMLSyntaxError, timeout=0.2)
    def parse_xml(self, xml_file):
        """Parse the XML file and populate the fields

        Retries 3 times if there is an XMLSyntaxError to avoid
        race condition flushing the file to disk

        Args:
            xml_file (str): Path to the XML file
        """
        return ElementTree.parse(xml_file).getroot()

    def _load_from_uniprot_xml(self, accession, xml_file_path=None) -> None:
        if xml_file_path is not None:
            xml_file = xml_file_path
        else:
            xml_file = fetch_uniprot_file(accession, "xml", unp_dir=self.unp_dir)
        doc = self.parse_xml(xml_file)

        accessions = list(map(str, doc[0].xpath(".//*[name()='accession']/text()")))
        if not accessions:
            raise ValueError(f"Invalid file for {accession}. No accessions found.")

        self._populate_fields(doc)

    def __init__(self, accession, unp_dir, fetch_annotation_score: bool = True):
        self.seq_isoforms: dict[str, Any] = {}
        self.isoforms: dict[str, Any] = {}
        self.features: dict[str, list] = {}
        self.comments: dict[str, list] = {}
        self.dbreferences: dict[str, list[str]] = {}
        self.evidences: list[dict] = []
        self.secondary_accessions: list[str] = []
        self.synonyms: list[str] = []
        self.accession = None
        self.annotation_score = None
        self.unp_dir = unp_dir
        self.ad_dbref_auto_acc = accession
        canonical = accession.split("-")[0] if "-" in accession else accession
        # Pass fetch_annotation_score so that when we save the pickle for the
        # first time we include the score — avoiding a separate network call on
        # every subsequent load.
        self._load_uniprot_data(canonical, fetch_annotation_score=fetch_annotation_score)
        # If the pickle was loaded but didn't yet contain annotation_score,
        # fetch it now and update the pickle.
        if fetch_annotation_score and self.annotation_score is None:
            self.annotation_score = get_annotation_score(canonical)
            self._save_to_pickle(self._pickle_path(canonical))

    def _populate_fields(self, doc):
        uniprot = doc[0]
        self.secondary_accessions = list(
            map(str, uniprot.xpath(".//*[name()='accession']/text()"))
        )
        self.accession = str(self.secondary_accessions[0])
        del self.secondary_accessions[0]

        about_sequence = uniprot.xpath(".//*[name()='sequence' and ./@length!='']")[0]
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

        if chainIndices != []:
            if (
                chainIndices[0].xpath(".//*[name()='begin']/@status") != []
                and chainIndices[0].xpath(".//*[name()='end']/@status") != []
            ):
                self.chain = (
                    int(chainIndices[0].xpath(".//*[name()='begin']/@position")[0]),
                    int(chainIndices[0].xpath(".//*[name()='end']/@position")[0]),
                )

        self.molName = [
            name.upper() for name in uniprot.xpath(".//*[name()='protein']//text()")
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

    def _extract_taxonomy(self, about_organism):
        self.organism = about_organism[0].xpath(
            ".//*[name()='name'][@type='scientific']/text()"
        ) + about_organism[0].xpath(".//*[name()='name'][@type='common']/text()")

        self.organism = list(map(str, self.organism))
        self.taxonomy = [
            str(tax.xpath("./@id")[0])
            for tax in about_organism[0].xpath(
                ".//*[name()='dbReference'][@type='NCBI Taxonomy']"
            )
        ]

    def _extract_features(self, about_features):
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
                else:  # some features have no @id in the XML (e.g. P26039)
                    id_value = None

                if feat.xpath(".//*[name()='original']/text()") == []:
                    original = ""
                    variation = ""
                else:
                    original = str(feat.xpath(".//*[name()='original']/text()")[0])
                    variation = str(feat.xpath(".//*[name()='variation']/text()")[0])

                if feat.xpath("./@evidence") == []:
                    evidence = None
                else:
                    evidence = str(feat.xpath("./@evidence")[0])

                location = feat.xpath(".//*[name()='location']")[0]

                if location.xpath(".//*[name()='position']") == []:
                    begin = location.xpath(".//*[name()='begin']/@position")[0]
                    end = location.xpath(".//*[name()='end']/@position")[0]
                else:
                    begin = location.xpath(".//*[name()='position']/@position")[0]
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
                        {"id": id_value, "evidence": evidence, "position": int(begin)}
                    )

    def _extract_comments(self, about_comments):
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

    def _extract_signal(self, signal_peptide):
        signal = None
        if signal_peptide:
            signal_peptide = signal_peptide[0]

            if (
                signal_peptide.xpath(".//*[name()='begin']/@status")
                and signal_peptide.xpath(".//*[name()='begin']/@status")[0] == "unknown"
            ):
                begin = None
            # if begin == end then they sometimes use position/position
            elif signal_peptide.xpath(".//*[name()='position']/@position"):
                begin = int(
                    signal_peptide.xpath(".//*[name()='position']/@position")[0]
                )
            else:
                begin = int(signal_peptide.xpath(".//*[name()='begin']/@position")[0])

            if (
                signal_peptide.xpath(".//*[name()='end']/@status")
                and signal_peptide.xpath(".//*[name()='end']/@status")[0] == "unknown"
            ):
                end = None
            elif signal_peptide.xpath(".//*[name()='position']/@position"):
                end = begin
            else:
                end = int(signal_peptide.xpath(".//*[name()='end']/@position")[0])

            signal = (begin, end)

        return signal

    def getAllIsoforms(self):
        if len(self.seq_isoforms) > 1:
            return self.seq_isoforms

        for iso, variants in list(self.isoforms.items()):
            # Don't follow isoforms with a different accession
            if iso.split("-")[0] == self.accession:
                self.seq_isoforms[iso] = self.applySpliceVariant(variants)

        return self.seq_isoforms

    def applySpliceVariant(self, spliceids):
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

    def getAllSpliceVariants(self):
        splices = []

        if self.hasSpliceVariants():
            splices.append(self.sequence)

            for var in self.features["splice variant"]:
                splices.append(self.applySpliceVariant(var["id"]))

        return splices

    def applySequenceVariant(self, variant):
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

    def getAllSequenceVariants(self):
        variants = []

        if self.hasVariants():
            variants.append(self.sequence)

            for var in self.features["sequence variant"]:
                variants.append(self.applySequenceVariant(var["id"]))

        return variants

    def isReviewed(self):
        return self.dataset == "Swiss-Prot"

    def getOrganismScientific(self):
        return self.organism if self.organism else None

    def getTaxId(self):
        return self.taxonomy if self.taxonomy else None

    def getUNPPrimaryName(self):
        return self.molName if self.molName else None

    def hasSecondaryName(self):
        return True if self.synonyms else False

    def getUNPSecondaryName(self):
        return self.synonyms if self.synonyms else []

    def hasShortName(self):
        return True if self.synonyms else False

    def getUNPShortName(self):
        return self.longName if self.longName else False

    def hasSpliceVariants(self):
        return (
            True
            if self.features and "splice variant" in list(self.features.keys())
            else False
        )

    def hasVariants(self):
        return (
            True
            if self.features and "sequence variant" in list(self.features.keys())
            else False
        )

    def getSpliceVariants(self):
        if self.hasSpliceVariants():
            return self.features["splice variant"]

        return []

    def getVariants(self):
        if self.hasVariants():
            return self.features["sequence variant"]

        return []

    def getSequenceOneLetter(self):
        return self.sequence

    def getAccession(self):
        return self.accession

    def getAccessionLong(self):
        return self.longName

    def hasSecondaryAccessions(self):
        return (
            True
            if self.secondary_accessions and self.secondary_accessions != []
            else False
        )

    def hasModifiedResidues(self):
        return (
            True
            if self.features and "modified residue" in list(self.features.keys())
            else False
        )

    def hasMutagenesisSite(self):
        return (
            True
            if self.features and "mutagenesis site" in list(self.features.keys())
            else False
        )

    def getFeatures(self, idx):
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

