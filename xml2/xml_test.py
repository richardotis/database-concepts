from pycalphad import Database, Model, calculate, variables as v
import xml_parser

dbf = Database('alzn_mey.xml')
#dbf = Database('alzn_mey.tdb')
mod = Model(dbf, ['AL', 'ZN'], 'LIQUID')
res = calculate(dbf, ['AL', 'ZN'], 'LIQUID', T=600, P=1e5)
print(res)

dbf.to_file('alzn_mey.xml', if_exists="overwrite")
