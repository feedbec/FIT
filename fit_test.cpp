#include<vector>
#include<iostream>

struct Node
{
    double glob_loc[3];
    std::vector<int> edges;
    Node(double x, double y, double z){glob_loc[0]=x;glob_loc[1]=y;glob_loc[2]=z;};

};
typedef std::vector<Node> grid;

struct Edge
{
int bnode;
int enode;
std::vector<int> faces;
Edge(int n1, int n2):bnode(n1),enode(n2){};
};
typedef std::vector<Edge> edge_grid;

struct Face
{
    std::vector<int> edges;
    int edge_dual; // if negative then this is the boundary face
    Face(){};
    Face(const std::vector<int>&a):edges(a){};
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
    Mesh m(4,4,4,0.1,0.1,0.1);
    // make nodes
    for(int i = 0; i<m.Nx*m.Ny*m.Nz; i++)
    {
        m.lin2vol(i,ix,iy,iz);
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
            m.common_edges(node_list, edge_list);
            m.fg.push_back(Face(edge_list));
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
            m.common_edges(node_list, edge_list);
            m.fg.push_back(Face(edge_list));
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
            m.common_edges(node_list, edge_list);
            m.fg.push_back(Face(edge_list));
            for (int j = 0; j < edge_list.size(); j++)
            {
                m.eg[edge_list[j]].faces.push_back(m.fg.size()-1);
            }
        }       
    }
    // find boundary and non-boundary edges
    std::vector<int> boundary_edges;
    std::vector<int> non_boundary_edges;

    for(int i = 0; i<m.eg.size();i++)
    {
        if(m.eg[i].faces.size() >= 4)
            non_boundary_edges.push_back(i);
        else
            boundary_edges.push_back(i);
    }
    // assemble matrix
    double M[m.eg.size()][m.eg.size()];
    for (int i = 0; i < m.eg.size(); i++)
    {
        for (int j = 0; j < m.eg.size(); j++)
        {
            M[i][j] = 0.0;
        }      
    }
    
    for(int i=0; i < m.eg.size();i++)
    {
        for(int j = 0; j < m.eg[i].faces.size(); j++)
        {
            for(int k = 0; k < m.fg[m.eg[i].faces[j]].edges.size(); k++)
                {
                    int edge_index = m.fg[m.eg[i].faces[j]].edges[k];
                    M[i][edge_index] += 1.0;
                }
        }
    }
    // remove boundary edges from matrix
    double M2[non_boundary_edges.size()][non_boundary_edges.size()];
    for(int i=0; i<non_boundary_edges.size(); i++)
        for(int j=0;j<non_boundary_edges.size();j++)
            M2[i][j] = M[non_boundary_edges[i]][non_boundary_edges[j]]; 

    
    for (int i = 0; i < non_boundary_edges.size(); i++)
    {
        for (int j = 0; j < non_boundary_edges.size(); j++)
        {
            std::cout<<M2[i][j]<<" ";
        }
        std::cout<<std::endl;      
    }
 return 0;
}