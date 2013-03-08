import os

from pyvmc.utils.torque import PBSFlags, get_array_ids

sample_file_contents = """garbage

#PBS -l walltime=24:00:00
#PBS -l nodes=1:ppn=4
#PBS -t 0-24

more garbage
"""

def test_pbs_flags(tmpdir):
    p = tmpdir.mkdir("sub").join("hello.txt")
    p.write(sample_file_contents)

    flags = PBSFlags(os.path.join(p.dirname, p.basename))

    assert flags.nodes == 1
    assert flags.ppn == 4
    assert flags.walltime == "24:00:00"

def test_array_ids():
    assert get_array_ids("20, 40, 60") == (20, 40, 60)
    assert get_array_ids("3") == (0, 1, 2)
    assert get_array_ids("1,4") == (1, 4)
    assert get_array_ids("1, 5-7, 8") == (1, 5, 6, 7, 8)
    assert get_array_ids("3-4") == (3, 4)
