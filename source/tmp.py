import diskhash

def f():
    s = diskhash.StructHash(fname='tmp.dht', keysize=10, structformat='llll', mode='rw')

    value=[1,2,3,4]
    print(s.insert('key4', *value))
    print(s.lookup('key'))
    print(s.lookup('key4'))
    print(s.size())

    print('abcd')

f()