import sketch
import fastx

p = sketch.Parameters()
msh1 = sketch.MinHash(p)
msh2 = sketch.MinHash(p)
for name, seq, qual in fastx.Fastx('genome1.fna'):
    msh1.update(seq)
for name, seq, qual in fastx.Fastx('genome2.fna'):
    msh2.update(seq)

jac = msh1.jaccard(msh2)

print("distance(1-JAC): ", 1.0-jac)
