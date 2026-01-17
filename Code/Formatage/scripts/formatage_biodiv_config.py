# formatage_biodiv_config.py

colonnes_import = {
    "GBIF": [
        "speciesKey", "species", "occurrenceStatus",
        "decimalLongitude", "decimalLatitude",
        "eventDate", "year", "month", "individualCount",
        "gridId", "grid_name", "taxonRank",
        "kingdom", "phylum", "class","order", "family", "genus"
    ],
    "INPN": [
        "cdRef", "nomScientifique", "occurrenceStatus",
        "x_centroid_4326", "y_centroid_4326",
        "dateObs", "year", "month", "nombreIndividus",
        "codeMaille10Km", "taxonRank",
        "kingdom", "phylum", "classe","ordre", "famille", "genre"
    ],
    "SILENE": [
        "cd_ref", "nom_valide","nom_vernaculaire",
        "x_centroid_4326", "y_centroid_4326",
        "date_debut", "nombre_min",
        "classe", "ordre", "famille","regne"
    ],
    "INAT": [
        "speciesKey", "species", "occurrenceStatus",
        "decimalLongitude", "decimalLatitude",
        "eventDate", "year", "month", "individualCount",
        "gridId", "grid_name", "taxonRank",
        "kingdom", "phylum", "class","order", "family", "genus"
    ],
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
    },
    "SILENE": {
        "cd_ref": "speciesID",
        "nom_valide": "nom_valide",
        "x_centroid_4326": "lon",
        "y_centroid_4326": "lat",
        "date_debut": "eventDate",
        "nombre_min": "nombreObs",
        "regne":"kingdom",
        "classe":"class",
        "ordre": "order",
        "famille": "family"
    },
    "INAT": {
        "speciesKey": "speciesID",
        "verbatimScientificName": "nom_valide",
        "decimalLongitude": "lon",
        "decimalLatitude": "lat",
        "dateIdentified": "eventDate",
        "kingdom":"kingdom",
        "phylum": "phylum",
        "class":"class",
        "order": "order",
        "family": "family",
        "genus": "genus",
        "occurrenceStatus":"occurrenceStatus",
        "taxonRank": "taxonRank"
    }
}
