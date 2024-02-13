#include<vector>
#include<iostream>
#include<cmath>
#include <fstream>
#include "Eigen/Sparse"

struct Node
{
    double glob_loc[3];
    std::vector<int> edges;
    double mu;
    Node(double x, double y, double z)
    {glob_loc[0]=x;glob_loc[1]=y;glob_loc[2]=z; mu = 1.0;};
    Node(double x, double y, double z, double mu0)
    {glob_loc[0]=x;glob_loc[1]=y;glob_loc[2]=z; mu = mu0;};

};
typedef std::vector<Node> grid;

struct Edge
{
int bnode;
int enode;
int sign;
int set_id;
double Val;
std::vector<int> faces;
Edge(int n1, int n2):bnode(n1),enode(n2)
{sign = 1; Val = 0.0; set_id = 0;};
Edge(){sign=1; Val = 0.0; set_id = 0;};
};
typedef std::vector<Edge> edge_grid;

struct Face
{
    std::vector<int> edges;
    double mu;
    Face(){mu = 1.0;};
    Face(const std::vector<int>&a, double mu0):edges(a){mu = mu0;};
};
typedef std::vector<Face> face_grid;

struct Cell
{
    std::vector<int> faces;
    std::vector<int> next_cells;
};

struct Mesh
{
    int Nx, Ny, Nz;
    double dx,dy,dz;
    grid g;
    edge_grid eg;
    face_grid fg;

    Mesh(int Nx0, int Ny0, int Nz0, double dx0, double dy0, double dz0):Nx(Nx0),
    Ny(Ny0),Nz(Nz0),dx(dx0),dy(dy0),dz(dz0){};
    void lin2vol(int j, int & ix, int & iy, int &iz);
    void vol2lin(int ix, int iy, int iz, int & j);
    void print_nodes();
    void print_edges();
    void print_faces();
    void common_edges(const std::vector<int>& node_list, std::vector<int>& edge_list);
};

void Mesh::print_faces()
{
    for (int i = 0; i < fg.size(); i++)
    {
        std::cout<< i <<" ";
        for (int j = 0; j < fg[i].edges.size(); j++)
        {
            std::cout<<fg[i].edges[j]<<" ";
        }
        std::cout<<std::endl;
    }
    
}
void Mesh::print_edges()
{
    for (int i = 0; i < eg.size(); i++)
    {
        std::cout<<i<<" nodes: "<<eg[i].bnode<<" "<<eg[i].enode<<" faces: ";
        for (int j = 0; j < eg[i].faces.size(); j++)
        {
            std::cout<<eg[i].faces[j]<<" ";
        }
        std::cout<<std::endl;
    }
    
}

void Mesh::print_nodes()
{
    for (int i = 0; i < g.size(); i++)
    {
        std::cout<<i<<" "<<g[i].glob_loc[0]<<" "
        <<g[i].glob_loc[1]<<" "<<g[i].glob_loc[2]<<" ";
        for(int j=0; j<g[i].edges.size(); j++)
            std::cout<<g[i].edges[j]<<" ";
        std::cout<<std::endl;
    }
    
}
void Mesh::common_edges(const std::vector<int>& node_list, std::vector<int>& edge_list)
{   
    edge_list.clear();
    std::vector<int> pool;
    for(int i = 0; i< g[node_list[0]].edges.size(); i++)
        pool.push_back(g[node_list[0]].edges[i]);
    
    for(int i=1; i<node_list.size();i++)
    {
        int pool_size = pool.size();
        for(int j = 0; j <g[node_list[i]].edges.size(); j++)
        {
            bool in_pool = false;
            for(int k = 0; k< pool_size; k++)
            {
                if (g[node_list[i]].edges[j] == pool[k])
                {
                    edge_list.push_back(pool[k]);
                    in_pool = true;
                    break;
                }
            }
            if(!in_pool)
            {
                pool.push_back(g[node_list[i]].edges[j]);
            }

        }
    }
}

void Mesh::lin2vol(int i, int & ix, int & iy, int &iz)
{
    ix= i/Nz/Ny;
    iy = (i - Nz*Ny*ix)/Nz;
    iz = i - Nz*Ny*ix - Nz*iy;
}
void Mesh::vol2lin(int ix, int iy, int iz, int & j)
{
    j = ix*Nz*Ny + iy*Nz + iz;
}
    
int main()
{
    int ix, iy, iz;
    int j;
    Mesh m(101,101,101,0.01,0.01,0.01);
    // make nodes
    for(int i = 0; i<m.Nx*m.Ny*m.Nz; i++)
    {   
        m.lin2vol(i,ix,iy,iz);
        if(ix*m.dx > 0.6)
            m.g.push_back(Node(m.dx*ix, m.dy*iy, m.dz*iz, 1.0));
        else
            m.g.push_back(Node(m.dx*ix, m.dy*iy, m.dz*iz));   
    }
    // make edges
    for(int i = 0; i<m.Nx*m.Ny*m.Nz; i++)
    {
        m.lin2vol(i,ix,iy,iz);
        if(ix+1<m.Nx){
            m.vol2lin(ix+1,iy,iz, j);
            m.eg.push_back(Edge(i,j));
            m.g[i].edges.push_back(m.eg.size()-1);
            m.g[j].edges.push_back(m.eg.size()-1);
        }
        if(iy+1<m.Ny){
            m.vol2lin(ix,iy+1,iz, j);
            m.eg.push_back(Edge(i,j));
            m.g[i].edges.push_back(m.eg.size()-1);
            m.g[j].edges.push_back(m.eg.size()-1);
        }
        if(iz+1<m.Nz){
            m.vol2lin(ix,iy,iz+1, j);
            m.eg.push_back(Edge(i,j));
            m.g[i].edges.push_back(m.eg.size()-1);
            m.g[j].edges.push_back(m.eg.size()-1);
        }
    }
    // make faces
    std::vector<int> node_list;
    std::vector<int> edge_list;
    int j1,j2,j3;
    for(int i= 0; i < m.Nx*m.Ny*m.Nz; i++)
    {
        node_list.clear();
        m.lin2vol(i,ix,iy,iz);
        if(ix+1<m.Nx && iy+1 < m.Ny){
            m.vol2lin(ix+1,iy,iz, j1);
            m.vol2lin(ix,iy+1,iz, j2);
            m.vol2lin(ix+1,iy+1,iz, j3);
            node_list.push_back(i);
            node_list.push_back(j1);
            node_list.push_back(j2);
            node_list.push_back(j3);
            double mu = 0.0;
            for(int k = 0; k < node_list.size(); k++)
                mu += m.g[node_list[k]].mu;
            mu = mu/node_list.size();
            m.common_edges(node_list, edge_list);
            m.fg.push_back(Face(edge_list, mu));
            for (int j = 0; j < edge_list.size(); j++)
            {
                m.eg[edge_list[j]].faces.push_back(m.fg.size()-1);
            }
        } 
        node_list.clear();
        if(ix+1<m.Nx && iz+1 < m.Nz){
            m.vol2lin(ix+1,iy,iz, j1);
            m.vol2lin(ix,iy,iz+1, j2);
            m.vol2lin(ix+1,iy,iz+1, j3);
            node_list.push_back(i);
            node_list.push_back(j1);
            node_list.push_back(j2);
            node_list.push_back(j3);
            double mu = 0.0;
            for(int k = 0; k < node_list.size(); k++)
                mu += m.g[node_list[k]].mu;
            mu = mu/node_list.size();
            m.common_edges(node_list, edge_list);
            m.fg.push_back(Face(edge_list, mu));
            for (int j = 0; j < edge_list.size(); j++)
            {
                m.eg[edge_list[j]].faces.push_back(m.fg.size()-1);
            }
        }       
        node_list.clear();
        if(iy+1<m.Ny && iz+1 < m.Nz){
            m.vol2lin(ix,iy+1,iz, j1);
            m.vol2lin(ix,iy,iz+1, j2);
            m.vol2lin(ix,iy+1,iz+1, j3);
            node_list.push_back(i);
            node_list.push_back(j1);
            node_list.push_back(j2);
            node_list.push_back(j3);
            double mu = 0.0;
            for(int k = 0; k < node_list.size(); k++)
                mu += m.g[node_list[k]].mu;
            mu = mu/node_list.size();
            m.common_edges(node_list, edge_list);
            m.fg.push_back(Face(edge_list, mu));
            for (int j = 0; j < edge_list.size(); j++)
            {
                m.eg[edge_list[j]].faces.push_back(m.fg.size()-1);
            }
        }       
    }
    // find boundary and non-boundary edges
    std::vector<int> boundary_edges;
    std::vector<int> non_boundary_edges;

    for(int i = 0; i<m.eg.size(); i++)
    {
        if(m.eg[i].faces.size() >= 4){
            non_boundary_edges.push_back(i);
            m.eg[i].set_id=non_boundary_edges.size()-1;// set the number of the edge in non_boundary array
        }
        else
            boundary_edges.push_back(i);
    }

    // assemble matrix
    std::cout<<"start matrix assembly..."<<std::endl;

    typedef Eigen::SparseMatrix<double> SpMat; // declares a column-major sparse matrix type of double
    typedef Eigen::Triplet<double> T;
    std::vector<T> coefficients;

    int next_node;
    std::vector<int> tmp_edge_array;
// run over non_boundary edges only
    for(int i_nd=0; i_nd < non_boundary_edges.size();i_nd++)
    {
        int i = non_boundary_edges[i_nd];
        next_node = m.eg[i].enode;
       // m.eg[i].sign = 1;
        coefficients.push_back(T(i_nd,i_nd,m.eg[i].faces.size()));
        for(int j = 0; j < m.eg[i].faces.size(); j++)
        {
            tmp_edge_array.clear();
            for(int k = 0; k < m.fg[m.eg[i].faces[j]].edges.size(); k++)
                {
                    int edge_index = m.fg[m.eg[i].faces[j]].edges[k];
                    if(edge_index != i)
                        tmp_edge_array.push_back(edge_index);
            }
            // determine edge orientation in the current face
            while (tmp_edge_array.size()>0)
            {
                for(int k = 0; k < tmp_edge_array.size(); k++)
                {
                    if(m.eg[tmp_edge_array[k]].bnode == next_node || m.eg[tmp_edge_array[k]].enode == next_node)
                    {
                        if(m.eg[tmp_edge_array[k]].bnode == next_node)
                        {
                            m.eg[tmp_edge_array[k]].sign = 1;
                            next_node = m.eg[tmp_edge_array[k]].enode;
                        }
                        else
                        {
                            m.eg[tmp_edge_array[k]].sign = -1;
                            next_node = m.eg[tmp_edge_array[k]].bnode;
                        }
                        tmp_edge_array.erase(tmp_edge_array.begin() + k); 
                        break;
                    }
                }
            }

            for(int k = 0; k < m.fg[m.eg[i].faces[j]].edges.size(); k++)
                {
                    int edge_index = m.fg[m.eg[i].faces[j]].edges[k];
                    if(edge_index != i)
                    {     
                        double coeff_value = m.eg[edge_index].sign/m.fg[m.eg[i].faces[j]].mu;
                        coefficients.push_back(T(i_nd,m.eg[edge_index].set_id, coeff_value));
                    }
                }
        }
    }
    SpMat Anb(non_boundary_edges.size(),non_boundary_edges.size());
    Anb.setFromTriplets(coefficients.begin(), coefficients.end());

    //set load vector
    Eigen::VectorXd b(non_boundary_edges.size()); 
    std::cout << "start loads..." << std::endl;
    for (int i = 0; i < non_boundary_edges.size(); i++)
    {
        double x0 = m.g[m.eg[non_boundary_edges[i]].bnode].glob_loc[0];
        double y0 = m.g[m.eg[non_boundary_edges[i]].bnode].glob_loc[1];
        double z0 = m.g[m.eg[non_boundary_edges[i]].bnode].glob_loc[2];
        double x1 = m.g[m.eg[non_boundary_edges[i]].enode].glob_loc[0];
        double y1 = m.g[m.eg[non_boundary_edges[i]].enode].glob_loc[1];
        double z1 = m.g[m.eg[non_boundary_edges[i]].enode].glob_loc[2];
        double pos = std::sqrt((x0-0.5)*(x0-0.5) + (y0-0.5)*(y0-0.5) );

        if(pos <0.005 && std::fabs(z0-z1) > 0.005)
        {
            if(z1>z0)
                b(i) = 1.0;
            else
                b(i) = -1.0;

            std::cout << i <<" "<<b(i) << std::endl;
        }        
        else
            b(i) = 0.0;
    }

 // solving
    std::cout << "start BiCGSTAB ..." << std::endl;
 
    Eigen::BiCGSTAB<SpMat> solver;
    solver.compute(Anb);

    Eigen::VectorXd x = solver.solve(b);
    std::cout << "#iterations:     " << solver.iterations() << std::endl;
    std::cout << "estimated error: " << solver.error()      << std::endl;
    for (int i = 0; i < non_boundary_edges.size(); i++)
    {
        m.eg[non_boundary_edges[i]].Val = x(i);
    }
    
    //write vector potential in file
    std::ofstream MyFile("Afile_no_mag.txt");

    for (int i = 0; i < m.g.size(); i++)
    {
        double Ax = 0.0;
        double Ay = 0.0;
        double Az = 0.0;

        for(int j = 0; j< m.g[i].edges.size();j++)
        {
            int edge_index = m.g[i].edges[j];
            double x0 = m.g[m.eg[edge_index].bnode].glob_loc[0];
            double y0 = m.g[m.eg[edge_index].bnode].glob_loc[1];
            double z0 = m.g[m.eg[edge_index].bnode].glob_loc[2];
            double x1 = m.g[m.eg[edge_index].enode].glob_loc[0];
            double y1 = m.g[m.eg[edge_index].enode].glob_loc[1];
            double z1 = m.g[m.eg[edge_index].enode].glob_loc[2];
            double L = std::sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0)); 
            Ax += std::abs(m.eg[edge_index].Val)*(x1-x0)/L;
            Ay += std::abs(m.eg[edge_index].Val)*(y1-y0)/L;
            Az += std::abs(m.eg[edge_index].Val)*(z1-z0)/L;
        }
        MyFile<<m.g[i].glob_loc[0]<<","<<m.g[i].glob_loc[1]<<","<<m.g[i].glob_loc[2]
        <<","<< Ax <<","<< Ay <<","<< Az <<std::endl; 
    }
    // Close the file
    MyFile.close();
    std::cout<<"Done!"<<std::endl;
    return 0;
}