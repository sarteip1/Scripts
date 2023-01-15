from pymongo import MongoClient

# source database
src_client = MongoClient("mongodb://ci1.fritz.box:27017/?readPreference=primary&directConnection=true&ssl=false")
src_db = src_client["DEFAULT_DB"]

# destination database
dst_client = MongoClient("mongodb://localhost:27017/?readPreference=primary&directConnection=true&ssl=false")
dst_db = dst_client["DEFAULT_DB"]

src_collections = src_db.list_collection_names()

for coll in src_collections:
    documents = src_db[coll].find()
    dst_db[coll].insert_many(documents)

src_client.close()
dst_client.close()
