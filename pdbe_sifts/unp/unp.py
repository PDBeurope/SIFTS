#!/usr/bin/env python3

"""
Filters information from UniProt XML and stores it in a dictionary.
This allows for faster access to a small subset of required information
from a large XML file.
"""

import os
import pickle
import random
import requests
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import funcy
from lxml import etree as ElementTree
from lxml.etree import XMLSyntaxError

from pdbe_sifts.base.log import logger
# from opensifts.base.pdbe_path import get_uniprot_cache_dir
from pdbe_sifts.base.utils import fetch_uniprot_file
# from release.config import Config

COLORS = {
    "white": 0,
    "red": 31,
    "green": 32,
    "yellow": 33,
    "blue": 34,
}

# conf = Config()
UNP_URL_SCORE = "https://rest.uniprot.org/uniprotkb/{uniprot_id}?fields=accession,annotation_score"

def get_annotation_score(uniprot_id: str) -> int | None:
    """
    Retrieve the annotation score of a UniProt entry.
    
    Args:
        uniprot_id (str): UniProt accession (e.g., "P05067")
    
    Returns:
        int | None: annotationScore if available, otherwise None
    """
    url = UNP_URL_SCORE.format(uniprot_id=uniprot_id)
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        score = data.get("annotationScore")
        if score:
            return score
        else:
            return 0
    else:
        print(f"Error {response.status_code} for {uniprot_id}")
        return None

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
    dbentry_id: str
    dbreferences: Mapping[str, list[str]] = {}
    features: Mapping[str, str] = {}
    dataset: str
    evidences: list[str] = []
    secondary_accessions: list[str] = []
    annotation_score: int

    def _load_uniprot_data(self, accession):
        """Load the data for the accession from the uniprot XML file

        Will try to load from a pickle file if it exists. If not, will fetch from uniprot and cache it.
        if env var ORC_NO_CACHE_ALL is set to True, will not cache the data nor load from cache.

        Args:
            accession (str): The accession to load

        """
        # cache = bool(os.getenv("ORC_NO_CACHE_ALL", True))

        # pickle_file_path = os.path.join(
        #     get_uniprot_cache_dir(accession), f"{accession}.pkl"
        # )

        # loaded = False
        # if cache and Path(pickle_file_path).exists():
        #     loaded = self._load_from_pickle(accession, pickle_file_path)
        #     if loaded:
        #         return

        self._load_from_uniprot_api(accession)

        # if cache:
        #     Path(pickle_file_path).parent.mkdir(parents=True, exist_ok=True)
        #     with open(pickle_file_path, "wb") as f:
        #         pickle.dump(self.__dict__, f, pickle.HIGHEST_PROTOCOL)

    # def _load_from_pickle(self, accession: str, pickle_file_path: str) -> bool:
    #     with open(pickle_file_path, "rb") as f:
    #         try:
    #             self.__dict__ = pickle.load(f)
    #             logger.info(f"Loaded UNP pickle from {pickle_file_path}")
    #             return True

    #         except Exception:
    #             logger.warning(
    #                 f"Could not load pickle file for {accession}. Will delete it and fetch from uniprot"
    #             )
    #             os.unlink(pickle_file_path)
    #     return False

    @funcy.retry(3, errors=XMLSyntaxError, timeout=0.2)
    def parse_xml(self, xml_file):
        """Parse the XML file and populate the fields

        Retries 3 times if there is an XMLSyntaxError to avoid
        race condition flushing the file to disk

        Args:
            xml_file (str): Path to the XML file
        """
        return ElementTree.parse(xml_file).getroot()

    def _load_from_uniprot_api(self, accession) -> dict:
        xml_file = fetch_uniprot_file(accession, "xml", unp_dir=self.unp_dir)
        doc = self.parse_xml(xml_file)

        accessions = list(map(str, doc[0].xpath(".//*[name()='accession']/text()")))
        if not accessions:
            raise ValueError(f"Invalid file for {accession}. No accessions found.")

        self.accession = accession
        accessions = list(map(str, doc[0].xpath(".//*[name()='accession']/text()")))
        if not accessions:
            raise ValueError(f"Invalid file for {accession}. No accessions found.")

        self.accession = accession
        self._populate_fields(doc)

    def __init__(self, accession, unp_dir):
        self.seq_isoforms: Mapping[str, Any] = {}
        self.isoforms: Mapping[str, Any] = {}
        self.features: Mapping[str, str] = {}
        self.accession = None
        self.unp_dir = unp_dir
        self.ad_dbref_auto_acc = accession
        if "-" in accession:
            accession = accession.split("-")[0]
        self._load_uniprot_data(accession)
        self.annotation_score = get_annotation_score(accession)

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
                else:  # check out P26039
                    id_value = str(random.randint(0, 1000))

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

            # splices.append(
            #     self.applySpliceVariant(
            #         [x["id"] for x in self.features["splice variant"]]
            #     )
            # )

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

