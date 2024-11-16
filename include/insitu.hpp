#ifndef INSITU_HPP
#define INSITU_HPP

#include "ascent.hpp"
#include "conduit_blueprint.hpp"


struct InSituVis {
    conduit::Node mesh;
    conduit::Node actions;

    ascent::Ascent a;

    InSituVis(size_t nx, size_t ny, double dx, double dy) {
        conduit::Node ascent_options;
        ascent_options["runtime/type"] = "ascent";
        ascent_options["runtime/vtkm/backend"] = "openmp";
        a.open(ascent_options);

        mesh["coordsets/coords/type"] = "uniform";
        mesh["coordsets/coords/dims/i"] = 2*nx;
        mesh["coordsets/coords/dims/j"] = ny;

        mesh["coordsets/coords/spacing/dx"] = dx;
        mesh["coordsets/coords/spacing/dy"] = dy;

        mesh["topologies/topo/type"] = "uniform";
        mesh["topologies/topo/coordset"] = "coords";

       
        std::cout << mesh.to_yaml() << std::endl;

        conduit::Node verify_info;
        if(!conduit::blueprint::mesh::verify(mesh, verify_info))
        {
            std::cout << "Mesh Verify failed!" << std::endl;
            std::cout << verify_info.to_yaml() << std::endl;
        }
        else
        {
            std::cout << "Mesh verify success!" << std::endl;
        }

        conduit::Node &add_act = actions.append();
        add_act["action"] = "add_scenes";

        conduit::Node& scenes = add_act["scenes"];
        scenes["s1/plots/p1/type"] = "pseudocolor";
        scenes["s1/plots/p1/field"] = "density";
        scenes["s1/image_prefix"] = "/home/samrusso/cambridge/cfd/assignment1/cfd/output/";
    }

    void write(auto&& rho, auto&& v1, auto&& v2, auto&& p, size_t index) {
        conduit::Node& scenes = actions.child(0); 
        mesh["fields/density/association"] = "vertex";
        mesh["fields/density/topology"] = "topo";

        mesh["fields/v1/association"] = "vertex";
        mesh["fields/v1/topology"] = "topo";

        mesh["fields/v2/association"] = "vertex";
        mesh["fields/v2/topology"] = "topo";

        mesh["fields/pressure/association"] = "vertex";
        mesh["fields/pressure/topology"] = "topo";

        mesh["fields/density/values"].set_external(rho);
        mesh["fields/v1/values"].set_external(v1);
        mesh["fields/v2/values"].set_external(v2);
        mesh["fields/pressure/values"].set_external(p);

        mesh["state/cycle"] = index;

        std::ostringstream oss;
        oss << "out" << index;
        scenes["s1/image_name"] = oss.str();

        // publish mesh to ascent
        a.publish(mesh);

        // execute the actions
        a.execute(actions);

        // std::cout << actions.to_yaml() << std::endl;

    }

    void close() {
        a.close();
    }
};






#endif // INSITU_HPP