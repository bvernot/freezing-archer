import sys, random
import time
start_time = time.time()
check_time = start_time

a = {}

dim = 1000

print "setting up %d arrays" % (dim*dim)

for i in range(dim):
    a[i] = {}
    for j in range(dim):
        a[i][j] = []
        pass
    pass


t = time.time()
print t-start_time, t-check_time
check_time = t



print "populating arrays"

for it in range(1000):
    for i in range(1000000):
        a[random.randrange(0,1000)][random.randrange(0,1000)].append(i)
        pass
    t = time.time()
    print "iteration %d | size=%d | total_t=%f | round_t=%f | len(a[0][0])=%s" % (it, sys.getsizeof(a), t-start_time, t-check_time, len(a[0][0]))
    check_time = t
    pass
