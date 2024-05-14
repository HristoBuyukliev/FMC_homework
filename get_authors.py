import pandas as pd
import numpy as np
from pydantic import BaseModel
from typing import List, Optional, Set, Dict
from Bio import Entrez
from datetime import datetime
from tqdm import tqdm
import unicodedata
import json
from collections import defaultdict


def normalize_name(name):
    # change e.g. Josien Régis to Josien Regis
    # or Tomažič Aleš to Tomazic Ales
    normalized_name = unicodedata.normalize("NFKD", name).encode("ASCII", "ignore").decode("utf-8")
    return normalized_name.lower()


def format_name_to_pubmed(person_profile):
    """
    convert the person_profiles.json style to pubmed style
    """
    if person_profile["lastName"] is None:
        return ""
    name = person_profile["lastName"].split(",")[0]
    name = name.split(" ")[0] + " " + name.split(" ")[-1]
    return normalize_name(name)


class Author(BaseModel):
    first_name: str
    last_name: str
    name: str
    affiliations: Set[str]
    articles: Set[str]
    match: Optional[bool]

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)


def generate_date_pairs(start_year=1960):
    """
    returns pairs of type ('1963/01/01', '1963/02/01')
    """
    current_year = datetime.now().year
    pairs = []
    for year in range(start_year, current_year + 1):
        for month in range(1, 12):
            pairs.append((f"{year}/{month:02d}/01", f"{year}/{month+1:02d}/01"))
        pairs.append((f"{year}/12/01", f"{year+1}/01/01"))
    return pairs


def search_pubmed(query, max_results=1_000_000, retmax=1000):
    # we do it per month, since pubmed has a limit of 9999 per query, even with pagination
    minmax_date_pairs = generate_date_pairs()
    ids = []
    Entrez.email = "hristo.buyukliev@gmail.com"
    for min_date, max_date in minmax_date_pairs:
        date_ids = []
        while len(ids) < max_results:
            handle = Entrez.esearch(
                db="pubmed", term=query, retmax=retmax, retstart=len(date_ids), mindate=min_date, maxdate=max_date
            )
            record = Entrez.read(handle)
            handle.close()
            new_ids = record["IdList"]
            date_ids += new_ids
            if len(new_ids) < retmax:
                break
        print(f"Between {min_date} and {max_date}, we found {len(date_ids)} articles")
        ids += date_ids
    return ids


def get_author_name(author: str):
    # human:
    first_name = author.get("ForeName", "").split(" ")[0]  # take care of the Alan B. Watts situation
    last_name = author.get("LastName", "")
    if first_name and last_name:
        return first_name, last_name
    # collective:
    collective_name = author.get("CollectiveName", None)
    if collective_name:
        return collective_name, ""
    raise Exception(f"Cant extract author name from author: {author}")


def process_record(record):
    article_info = {}
    try:
        article_info["title"] = record["MedlineCitation"]["Article"]["ArticleTitle"]
    except KeyError:
        article_info["title"] = "No title available"
    article_info["pmid"] = str(record["MedlineCitation"]["PMID"].strip())
    authors = record["MedlineCitation"]["Article"].get("AuthorList", [])
    author_list = []
    for author in authors:
        first_name, last_name = get_author_name(author)
        author_affiliations = author.get("AffiliationInfo", [])
        author_affiliations = set([aa.get("Affiliation", "") for aa in author_affiliations])
        author_info = Author(
            first_name=first_name,
            last_name=last_name,
            name=f"{first_name} {last_name}",
            affiliations=author_affiliations,
            articles={article_info["pmid"]},
        )
        author_list.append(author_info)
    article_info["authors"] = author_list
    pubmed_data = record.get("PubmedData", {})
    history = pubmed_data.get("History", [])
    for entry in history:
        if entry.attributes["PubStatus"] == "pubmed":
            pub_date = entry
            break
    else:
        pub_date = None
    if pub_date:
        publication_date = "{}/{}/{}".format(pub_date["Month"], pub_date["Day"], pub_date["Year"])
    else:
        publication_date = "Unknown"
    article_info["publication_date"] = publication_date
    return article_info


def fetch_article_info(pmids, batch_size=1000):
    articles = []
    article_batches = [pmids[i : i + batch_size] for i in range(0, len(pmids), batch_size)]
    for batch in tqdm(article_batches):
        handle = Entrez.efetch(db="pubmed", id=",".join(batch), retmode="xml")
        try:
            records = Entrez.read(handle)
        except:
            print("problem in batch, skipping...")
            continue
        articles += records["PubmedArticle"] + records["PubmedBookArticle"]
    handle.close()
    processed_articles = []
    for article in articles:
        try:
            processed_articles.append(process_record(article))
        except:
            pass
    return processed_articles


query = "Ulcerative Colitis"
pmids = search_pubmed(query)
article_infos = fetch_article_info(pmids)

# get all authors:
authors = [author for article_info in article_infos for author in article_info["authors"]]
print(f"We found {len(authors)} authors.")


# deduplicate
def authored_by(article_info, author_name):
    for author in article_info["authors"]:
        if author_name == author["name"]:
            return True
    return False


deduped_authors = set(authors)
print(f"After deduplication, we are left with {len(deduped_authors)} authors")

# make sure that each author has ALL of his article ids; and all of his affiliations
# we'll do it linearly over the article_infos
author_dict: Dict[str, Author] = {author.name: author for author in deduped_authors}
for article in article_infos:
    for author in article["authors"]:
        author_dict[author.name].articles.add(article["pmid"])
        author_dict[author.name].affiliations.update(author.affiliations)
deduped_authors = list(author_dict.values())

# Finally, see how many of these appear in the provided dataset:
with open("person_profiles.json", "r") as rfile:
    person_profiles = json.load(rfile)
person_profile_names = [format_name_to_pubmed(pp) for pp in person_profiles]


matches = 0
for author in deduped_authors:
    cleaned_name = normalize_name(author.name)
    found_match = False
    if cleaned_name in person_profile_names:
        found_match = True
    author.match = found_match
    if found_match:
        matches += 1
        continue
print(f"We found a good match for {matches} out of {len(deduped_authors)} authors")

resulting_data = pd.DataFrame(
    {
        "fullName": [author.name for author in deduped_authors],
        "firstName": [author.first_name for author in deduped_authors],
        "lastName": [author.last_name for author in deduped_authors],
        "affiliations": [list(author.affiliations) for author in deduped_authors],
        "article_ids": [list(author.articles) for author in deduped_authors],
        "foundMatch": [author.match for author in deduped_authors],
    }
)

resulting_data.to_csv("pubmed_scrape.csv", index=False)
