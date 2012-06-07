import sys

# Set path to library here.
pathToLib = '/Users/paul/work/packages/irtk/build-64/lib'

sys.path.append(pathToLib)


import irtk

attr = irtk.irtkImageAttributes()

attr._x = 5
attr._y = 3
attr._z = 4

for i in range(20):
    print attr.IndexToLattice(i)

print
l = attr.IndexToLattice(26)

print l
print len(l)

im = irtk.irtkGreyImage()
im.Initialize(attr)
im.Print()

print im.ImageToWorld(0, 0, 0)

print im.WorldToImage(0, 0, 0)

print

print "image orientation"
print im.Orientation()


print



tr = irtk.irtkRigidTransformation()
tr.PutTranslationX(2)
tr.PutTranslationY(-1)
tr.PutRotationZ(90)

tr.Print()

pt = irtk.irtkPoint()
pt._x = 3

print
print "image to world on point"
print pt
im.ImageToWorld(pt)
print pt
print

print
print "world to image on point"
print pt
im.WorldToImage(pt)
print pt
print

print tr.Transform(3, 0, 0)

print
print pt
tr.Transform_point(pt)
print pt


af = irtk.irtkAffineTransformation()

af.Print()
af.PutScaleX(200)
af.Print()

af.Transform_point(pt)
print pt

mat = af.GetMatrix()
print "Matrix:"
print mat
print


print " print dir(irtk.irtkTransformation) "
print dir(irtk.irtkTransformation)
print "\n\n"

print "dir(irtk.irtkFreeFormTransformation)"
print dir(irtk.irtkFreeFormTransformation)
print "\n\n"

print "dir(irtk.irtkFreeFormTransformation3D)"
print dir(irtk.irtkFreeFormTransformation3D)
print "\n\n"

print "dir(irtk.irtkBSplineFreeFormTransformation3D)"
print dir(irtk.irtkBSplineFreeFormTransformation3D)
print "\n\n"

print 
print dir(irtk.irtkBSplineFreeFormTransformation3D.Print)

bffd3 = irtk.irtkBSplineFreeFormTransformation3D()
print
print "BFFD3D"
bffd3.Print()

print bffd3


mffd = irtk.irtkMultiLevelFreeFormTransformation()
mffd.Print()

print mffd.NumberOfLevels()

print dir(mffd)


# Following still breaks, cannot append a local transformation.
# mffd.PushLocalTransformation(bffd3)


