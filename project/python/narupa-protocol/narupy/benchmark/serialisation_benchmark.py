from narupa.protocol.benchmark.benchmark_pb2 import RawFrame, RawBytes
import time
import statistics
import random
import ctypes

for count in [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000]:

    times = []

    iters = int(10000000 / count)

    for i in range(iters):
        data = [random.uniform(0, 1) for _ in range(count)]
        frame = RawBytes()
        buf = (ctypes.c_float * count)()

        t0 = time.time()

        buf[:] = data
        frame.bytes = bytes(buf)
        output = frame.SerializeToString()

        t1 = time.time()

        times.append((t1 - t0) * 1000000)


    print("%f %f %f %i" % (statistics.median(times), statistics.mean(times), statistics.stdev(times), count))

