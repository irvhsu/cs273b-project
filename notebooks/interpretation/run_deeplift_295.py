import subprocess

num = 200
procs = 5
out = lambda x: "deeplift_out/deep_test{}.csv".format(x)
frags = range(0, num, num / procs) + [num]

tasks = []
for start, end in zip(frags[:-1], frags[1:]):
    cmd = "python get_deeplift_295.py {} {} {}".format(start, end, out(start))
    print cmd
    tasks += [subprocess.Popen(cmd, shell = True)]

for task in tasks: task.wait()
