#!/usr/bin/env python

"""
Run once for initial setup of db namespace.

"""

from BioSQL import BioSeqDatabase

# Set up the server connection details
server = BioSeqDatabase.open_database(driver = "MySQLdb", user = "root", passwd = "mmgdatabase", host = "localhost", db = "mmg_test")

db = server.new_database("mmg", description = "Testing MMG")

server.commit()