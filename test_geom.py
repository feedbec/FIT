# ellpsoid with holes
import pygmsh

with pygmsh.occ.Geometry() as geom:
    geom.characteristic_length_max = 0.1
    ellipsoid = geom.add_ellipsoid([0.0, 0.0, 0.0], [1.0, 0.7, 0.5])

    cylinders = [
        geom.add_cylinder([-1.0, 0.0, 0.0], [2.0, 0.0, 0.0], 0.3),
        geom.add_cylinder([0.0, -1.0, 0.0], [0.0, 2.0, 0.0], 0.3),
        geom.add_cylinder([0.0, 0.0, -1.0], [0.0, 0.0, 2.0], 0.3),
    ]
    geom.boolean_difference(ellipsoid, geom.boolean_union(cylinders))

    mesh = geom.generate_mesh()
    pygmsh.write("test1.msh")
    
