"""
An example in which 3 functions of x and y  are displayed with a surf plot,
while the z scaling is kept constant, to allow comparison between them.

The important aspect of this example is that the 3 functions should not
be displayed on top of each other, but side by side. For this we use the
extent keyword argument.

In addition, the relative scale between the different plots is important.
This is why we also use the `warp_scale` keyword argument, to have the same
scale on all plots.

Finally, we have to adjust the data bounds: as we want the "horizon" of
the wigner function in the middle of our extents, we put this to zero.

We add a set of axes and outlines to the plot. We have to play we extents
and ranges in order to make them fit with the data.
"""
# Author: Gael Varoquaux <gael.varoquaux@normalesup.org>
# Copyright (c) 2007, Enthought, Inc.
# License: BSD Style.

"""
from apptools.scripting.api import Recorder, set_recorder
# Create a recorder.
r = Recorder()
# Set the global recorder so the decorator works.
set_recorder(r)
r.register(p)
r.recording = True
"""


import WrightTools as wt
import numpy
import numpy as np
from mayavi import mlab
from mayavi.sources.builtin_surface import BuiltinSurface
from mayavi.modules.surface import Surface
from mayavi.filters.transform_data import TransformData

"""
data = wt.data.from_pickle(r'C:\Users\Loaner\Box Sync\Wright Shared\MX2\2015.02.08\Post-optimization data\Fluence Study\data.p')
data.level(0,'d2',-3)
data.zoom(3)
data.smooth([6,6,0])
data.convert('eV')
"""

"""
dims = np.array((23,41,41))
vol = np.array((data.d2.points.min(), data.d2.points.max(), data.w1.points.min(), data.w1.points.max(), data.w2.points.max(), data.w2.points.min()))
origin = vol[::2]
spacing = (vol[1::2]-origin)/(dims-1)
xmin, xmax, ymin, ymax, zmin, zmax = vol
x, y, z = np.ogrid[xmin:xmax:dims[0]*1j, ymin:ymax:dims[2]*1j, zmin:zmax:dims[1]*1j]
s = data.channels[0].values[18,:,:]


cat3_extent = (0, 41, 0, 41, 0, 60)
surf_proj = mlab.contour3d(s, colormap='spectral',
            extent=cat3_extent, vmin=data.w1.points.min(), vmax=data.w1.points.max())
mlab.outline(surf_proj, color=(.7, .7, .7), extent=cat3_extent)

#mlab.text(6, -2.5, '3 photons', z=-4, width=0.14)

#mlab.title('Multi-photons cats Wigner function')


dims = np.array((41,41,23))
vol = np.array((data.w2.points.min(), data.w2.points.max(), data.w1.points.min(), data.w1.points.max(), data.d2.points.max(), data.d2.points.min()))
origin = vol[::2]
spacing = (vol[1::2]-origin)/(dims-1)
xmin, xmax, ymin, ymax, zmin, zmax = vol
x, y, z = np.ogrid[xmin:xmax:dims[0]*1j, ymin:ymax:dims[1]*1j, zmin:zmax:dims[2]*1j]
s = data.channels[0].values[:,24,:].T


cat3_extent = (0, 41, 0, 41, 0, -60)
surf_proj = mlab.contour3d(s, colormap='spectral',
            extent=cat3_extent, vmin=data.w1.points.min(), vmax=data.w1.points.max())
mlab.outline(surf_proj, color=(.7, .7, .7), extent=cat3_extent)
"""
#mlab.text(6, -2.5, '3 photons', z=-4, width=0.14)


"""
dims = np.array((41,23,41))
vol = np.array((data.w2.points.min(), data.w2.points.max(), data.w1.points.min(), data.w1.points.max(), data.d2.points.max(), data.d2.points.min()))
origin = vol[::2]
spacing = (vol[1::2]-origin)/(dims-1)
xmin, xmax, ymin, ymax, zmin, zmax = vol
x, y, z = np.ogrid[xmin:xmax:dims[0]*1j, ymin:ymax:dims[1]*1j, zmin:zmax:dims[2]*1j]
s = data.channels[0].values[24,:,:]


cat3_extent = (0, 41, 41, 0, 0, -100)
surf_proj = mlab.contour3d(s, colormap='spectral',
            extent=cat3_extent, vmin=data.w1.points.min(), vmax=data.w1.points.max())
mlab.outline(surf_proj, color=(.7, .7, .7), extent=cat3_extent)

mlab.text(6, -2.5, '3 photons', z=-4, width=0.14)
"""

"""
dims = np.array((41,23,41))
vol = np.array((data.w2.points.min(), data.w2.points.max(), data.w1.points.min(), data.w1.points.max(), data.d2.points.max(), data.d2.points.min()))
origin = vol[::2]
spacing = (vol[1::2]-origin)/(dims-1)
xmin, xmax, ymin, ymax, zmin, zmax = vol
x, y, z = np.ogrid[xmin:xmax:dims[0]*1j, ymin:ymax:dims[1]*1j, zmin:zmax:dims[2]*1j]
s = data.channels[0].values[:,:,:].T
#mlab.figure(bgcolor=(1, 1, 1))
mlab.pipeline.volume(src, vmin=0.0999, vmax=0.10)
#mlab.contour3d(x, y, z, lab.pipeline.scalar_field(s), vmin=0.05, vmax=0.15, color=(1,0,0))
"""
"""
data = wt.data.from_pickle(r'C:\Users\Loaner\Box Sync\Wright Shared\MX2\2015.02.08\Post-optimization data\Fluence Study\ND2=0.2\data.p')
data.level(0,'d2',-3)
data.zoom(4)
data.smooth([8,8,0])
data.convert('eV')

dims = np.array((len(data.w2.points),len(data.w2.points),len(data.d2.points)))
vol = np.array((data.w2.points.min(), data.w2.points.max(), data.w1.points.min(), data.w1.points.max(), data.d2.points.max(), data.d2.points.min()))
origin = vol[::2]
spacing = (vol[1::2]-origin)/(dims-1)
xmin, xmax, ymin, ymax, zmin, zmax = vol
x, y, z = np.ogrid[xmin:xmax:dims[0]*1j, ymin:ymax:dims[1]*1j, zmin:zmax:dims[2]*1j]
s = data.channels[0].values[:,:,:]
#mlab.figure(bgcolor=(1, 1, 1))
src = mlab.pipeline.scalar_field(s)
#mlab.pipeline.volume(src, vmin=0.0999, vmax=0.10)
mlab.pipeline.iso_surface(src, contours=[s.min()+0.3*s.ptp(), ], opacity=0.3)
"""

data = wt.data.from_pickle(r'C:\Users\Loaner\Box Sync\Wright Shared\MX2\2015.02.08\Post-optimization data\Fluence Study\ND2=0.2\data.p')
data.level(0,'d2',-3)
data.zoom(3)
data.smooth([6,6,0])
data.convert('eV')



dims = np.array((len(data.w2.points),len(data.w2.points),len(data.d2.points)))
vol = np.array((data.w2.points.min(), data.w2.points.max(), data.w1.points.min(), data.w1.points.max(), data.d2.points.max(), data.d2.points.min()))
origin = vol[::2]
spacing = (vol[1::2]-origin)/(dims-1)
xmin, xmax, ymin, ymax, zmin, zmax = vol
x, y, z = np.ogrid[xmin:xmax:dims[0]*1j, ymin:ymax:dims[1]*1j, zmin:zmax:dims[2]*1j]
s = data.channels[0].values[:,:,:]
s1 = data.channels[0].values[:,:,:]
s2 = data.channels[0].values[:,:,:]
s3 = data.channels[0].values[:,:,:]
mlab.figure(bgcolor=(1, 1, 1))
src = mlab.pipeline.scalar_field(s, origin = (2.0, 0.0,100.0))
src1 = mlab.pipeline.scalar_field(s1)
src2 = mlab.pipeline.scalar_field(s2)
src3 = mlab.pipeline.scalar_field(s3)
#mlab.pipeline.volume(src, vmin=0.0999, vmax=0.10)
"""
mayavi.filters.image_change_information(src)
"""
"""
Instance(tvtk.ImageChangeInformation, args=(),
                      allow_none=False, record=True)
"""


planes1 = mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                            plane_orientation='x_axes',
                            slice_index=10,
                        )
                    
           
               
planes2 = mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                            plane_orientation='y_axes',
                            slice_index=10,
                        )

planes3 = mlab.pipeline.image_plane_widget(mlab.pipeline.scalar_field(s),
                            plane_orientation='z_axes',
                            slice_index=23,
                        )


#planes1.implicit_plane.widget.enabled = False
#planes2.implicit_plane.widget.enabled = False
#planes3.implicit_plane.widget.enabled = False
surf = mlab.pipeline.iso_surface(src, contours=[s.max()-0.4*s.ptp(), ],)
"""
#planes5 = mlab.pipeline.scalar_field(s, origin=(1.0,1.0,10.0))
planes4 = mlab.pipeline.scalar_cut_plane(mlab.pipeline.scalar_field(s),
                            plane_orientation='z_axes',
                        )

planes4.implicit_plane._set_origin((1.0, 1.0, 3.0))
planes4.implicit_plane.widget.enabled = False
"""
"""
lut = surf.module_manager.scalar_lut_manager.lut.table.to_array()


lut = rgba*255

planes1.module_manager.scalar_lut_manager.lut.table = lut

planes1 = mlab.pipeline.iso_surface(src, contours=[s.max()-0.4*s.ptp(), ],)

lut = planes1.module_manager.scalar_lut_manager.lut.table.to_array()


lut = rgba*255

planes1.module_manager.scalar_lut_manager.lut.table = lut


planes2.module_manager.scalar_lut_manager.lut.table = lut

planes2 = mlab.pipeline.iso_surface(src, contours=[s.max()-0.4*s.ptp(), ],)

lut = planes2.module_manager.scalar_lut_manager.lut.table.to_array()


lut = rgba*255

planes2.module_manager.scalar_lut_manager.lut.table = lut

planes3.module_manager.scalar_lut_manager.lut.table = lut

planes3 = mlab.pipeline.iso_surface(src, contours=[s.max()-0.4*s.ptp(), ],)

lut = planes3.module_manager.scalar_lut_manager.lut.table.to_array()


lut = rgba*255

planes3.module_manager.scalar_lut_manager.lut.table = lut
"""

"""

surf1 = mlab.pipeline.iso_surface(src1, contours=[s1.max()-0.4*s1.ptp(), ],)

lut1 = surf.module_manager.scalar_lut_manager.lut.table.to_array()


lut1 = rgba*255

surf1.module_manager.scalar_lut_manager.lut.table = lut1
"""

#mayavi.tools.pipeline.scalar_field(s)
#mayavi.tools.pipeline.scalar_field(s)
#mayavi.tools.pipeline.scalar_field(s)




#mlab.outline()

"""
n=3
fakear = np.zeros((len(data.w2.points)*n,len(data.w1.points)*n,len(data.d2.points)*n))
vects = np.array(([2,0,1],[2,1,2],[0,1,1], [1,1,1]))
lens = np.array([len(data.w2.points), len(data.w1.points), len(data.d2.points)])
for p in range(4):
    for i in range(len(data.w2.points)):
        for j in range(len(data.w1.points)):
            for k in range(len(data.d2.points)):
                #print lens*vects[p][0]
                #print lens*vects[p][1]
                #print lens*vects[p][2]
                fakear[i+(lens*vects[p][0])[0], j+(lens*vects[p][0])[1], k+(lens*vects[p][0])[2]]= s[i,j,k]
"""
"""
dims = np.array((41,41,23))
vol = np.array((data.d2.points.min(), data.d2.points.max(), data.w1.points.min(), data.w1.points.max(), data.w2.points.max(), data.w2.points.min()))
origin = vol[::2]
spacing = (vol[1::2]-origin)/(dims-1)
xmin, xmax, ymin, ymax, zmin, zmax = vol
x, y, z = np.ogrid[xmin:xmax:dims[0]*1j, ymin:ymax:dims[2]*1j, zmin:zmax:dims[1]*1j]
s = data.channels[0].values[18,:,:].T


cat3_extent = (0, 41, 0, 41, 0, 60)
surf_proj = mlab.contour3d(s, colormap='spectral',
            extent=cat3_extent, vmin=data.w1.points.min(), vmax=data.w1.points.max())
mlab.outline(surf_proj, color=(.7, .7, .7), extent=cat3_extent)

dims = np.array((41,41,23))
vol = np.array((data.d2.points.min(), data.d2.points.max(), data.w1.points.min(), data.w1.points.max(), data.w2.points.max(), data.w2.points.min()))
origin = vol[::2]
spacing = (vol[1::2]-origin)/(dims-1)
xmin, xmax, ymin, ymax, zmin, zmax = vol
x, y, z = np.ogrid[xmin:xmax:dims[0]*1j, ymin:ymax:dims[2]*1j, zmin:zmax:dims[1]*1j]
s = data.channels[0].values[18,:,:]


cat3_extent = (0, 41, 0, 41, 0, 60)
surf_proj = mlab.contour3d(s, colormap='spectral',
            extent=cat3_extent, vmin=data.w1.points.min(), vmax=data.w1.points.max())
mlab.outline(surf_proj, color=(.7, .7, .7), extent=cat3_extent)
"""

"""
src = mlab.pipeline.scalar_field(s)
mlab.pipeline.iso_surface(src, contours=[s.min()+0.1*s.ptp(), ], opacity=0.1)
mlab.pipeline.iso_surface(src, contours=[s.max()-0.1*s.ptp(), ],)
mlab.pipeline.image_plane_widget(src,
                            plane_orientation='z_axes',
                            slice_index=10,
                        )


"""
"""
def rotMat3D(axis, angle, tol=1e-12):
    #Return the rotation matrix for 3D rotation by angle `angle` degrees about an
    #arbitrary axis `axis`.
    
    t = np.radians(angle)
    x, y, z = axis
    R = (np.cos(t))*np.eye(3) +\
    (1-np.cos(t))*np.matrix(((x**2,x*y,x*z),(x*y,y**2,y*z),(z*x,z*y,z**2))) + \
    np.sin(t)*np.matrix(((0,-z,y),(z,0,-x),(-y,x,0)))
    R[np.abs(R)<tol]=0.0
    return R

# Main code

fig = mlab.figure()

engine = mlab.get_engine()

# Add a cylinder builtin source
cylinder_src = BuiltinSurface()
engine.add_source(cylinder_src)
cylinder_src.source = 'cylinder'
cylinder_src.data_source.center = np.array([ 0.,  0.,  0.])
cylinder_src.data_source.radius = 1.0
cylinder_src.data_source.capping = False
cylinder_src.data_source.resolution = 25

# Add transformation filter to rotate cylinder about an axis
transform_data_filter = TransformData()
engine.add_filter(transform_data_filter, cylinder_src)
Rt = np.eye(4)
Rt[0:3,0:3] = rotMat3D((1,0,0), 0) # in homogeneous coordinates
Rtl = list(Rt.flatten()) # transform the rotation matrix into a list

transform_data_filter.transform.matrix.__setstate__({'elements': Rtl})
transform_data_filter.widget.set_transform(transform_data_filter.transform)
transform_data_filter.filter.update()
transform_data_filter.widget.enabled = False   # disable the rotation control further.

# Add surface module to the cylinder source
cyl_surface = Surface()
engine.add_filter(cyl_surface, transform_data_filter)
# add color property
cyl_surface.actor.property.color = (1.0, 0.0, 0.0)

mlab.show()
"""
"""
mlab.draw()
mlab.view(40, 85)
"""
mlab.show()

#mlab.xlabel("X")
#mlab.show()
