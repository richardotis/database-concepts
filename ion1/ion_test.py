from pycalphad import Database, Model, calculate, variables as v
import ion_parser

dbf_ion = Database('test.ion')
mod = Model(dbf_ion, ['AL', 'ZN'], 'LIQUID')
res = calculate(dbf_ion, ['AL', 'ZN'], 'LIQUID', T=600, P=1e5)
print(res)