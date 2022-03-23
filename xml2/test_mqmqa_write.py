from pycalphad import Database, Model, calculate, variables as v
from pycalphad.tests.fixtures import select_database, load_database
import xml_parser

@select_database("Shishin_Fe-Sb-O-S_slag.dat")
def test_write_shishin_xml(load_database):
    """Test writing a DAT file with an MQMQA model validates correctly when writing"""
    dbf = load_database()
    # If the database fails to validate, it should raise
    dbf.to_string(fmt="xml", if_exists="overwrite")


if __name__ == "__main__":
    from importlib_resources import files
    from pycalphad.tests import databases
    dbf = Database(str(files(databases) / "Shishin_Fe-Sb-O-S_slag.dat"))
    mod = Model(dbf, ['FE', 'O', 'VA'], 'SLAG-LIQ')
    res = calculate(dbf, ['FE', 'O', 'VA'], 'SLAG-LIQ', T=600, P=1e5)
    print(res)
    dbf.to_file('Shishin_MQMQA.xml', if_exists="overwrite")

