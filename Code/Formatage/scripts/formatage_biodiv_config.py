# formatage_biodiv_config.py

# Colonnes à importer selon la source
colonnes_import = {
    "GBIF": [
        "speciesKey", "species", "occurrenceStatus",
        "decimalLongitude", "decimalLatitude",
        "eventDate", "year", "month", "individualCount",
        "gridId", "grid_name", "taxonRank",
        "kingdom", "phylum", "order", "family", "genus"
    ],
    "INPN": [
        "cdRef", "nomScientifique", "occurrenceStatus",
        "longitude", "latitude",
        "dateObs", "year", "month", "nombreIndividus",
        "codeMaille10Km", "taxonRank",
        "kingdom", "phylum", "ordre", "famille", "genre"
    ]
}

# Correspondance colonnes source → standard
col_map = {
    "GBIF": {
        "speciesKey": "speciesID",
        "species": "species",
        "decimalLongitude": "lon",
        "decimalLatitude": "lat",
        "eventDate": "eventDate",
        "year": "year",
        "month": "month",
        "individualCount": "nombreObs",
        "gridId": "gridId",
        "grid_name": "grid_name",
        "taxonRank": "taxonRank",
        "kingdom": "kingdom",
        "phylum": "phylum",
        "class":"class",
        "order": "order",
        "family": "family",
        "genus": "genus",
        "occurrenceStatus":"occurrenceStatus"
    },
    "INPN": {
        "cdRef": "speciesID",
        "nomScientifique": "species",
        "longitude": "lon",
        "latitude": "lat",
        "dateObs": "eventDate",
        "year": "year",
        "month": "month",
        "nombreIndividus": "nombreObs",
        "codeMaille10Km": "grid_name",
        "taxonRank": "taxonRank",
        "kingdom": "kingdom",
        "phylum": "phylum",
        "class":"class",
        "ordre": "order",
        "famille": "family",
        "genre": "genus"
    }
}
