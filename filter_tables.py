"""
    Filter the table for specific criteria
    Specifically, we want autism/control phenotypes, and SNV only. 
    Also get positions, chromosomes, genes. Other fields are generally
        irrelevant. No intergenic SNVs.
"""


import pymysql
import getpass


def db_filter(cur):
    """
        Filter table given cursor object to MySQL database
    """
    # fields to extract
    selected_fields = [
        "PrimaryPhenotype", 
        "Gene",
        "Transcript",
        "Chromosome", 
        "Position", 
        "Variant" 
    ]
    cur.execute(
        "SELECT {} FROM denovo_db".format(','.join(selected_fields))
        + " WHERE (PrimaryPhenotype='autism"
        + " OR PrimaryPhenotype='control')"
        + " AND Variant RLIKE '^[A-Z]>[A-Z]$';"
    )
    return cur.fetchall()


def main():
    """
        Main function
    """
    connection = pymysql.connect(
        host='bm185s-mysql.ucsd.edu',
        user='root',
        db='bimm185',
        passwd=getpass.getpass("Input password: ")
    )
    cursor = connection.cursor()
    records = db_filter(cursor)
