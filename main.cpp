#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/per_vertex_normals.h>
#include <igl/outer_element.h>
#include <igl/per_face_normals.h>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <igl/edge_lengths.h>


using namespace std;
using namespace Eigen;

namespace mesh
{
	template <typename DerivedV, typename DerivedF, typename DeriveddblA>
	IGL_INLINE void doublearea(
	const Eigen::MatrixBase<DerivedV> & V,
	const Eigen::MatrixBase<DerivedF> & F,
	Eigen::PlainObjectBase<DeriveddblA> & dblA)
	{
		const int dim = V.cols();
		const size_t m = F.rows();

		static int debugflag = 1;
		// Compute edge lengths
		//Eigen::Matrix<typename DerivedV::Scalar, Eigen::Dynamic, 3> l;
		//igl::edge_lengths(V,F,l);

		//igl::doublearea(l,0.0,dblA);

		// Projected area helper
		const auto & proj_doublearea =
		[&V,&F](const int x, const int y, const int f)
		->typename DerivedV::Scalar
		{
			auto rx = V(F(f,0),x)-V(F(f,2),x);
			auto sx = V(F(f,1),x)-V(F(f,2),x);
			auto ry = V(F(f,0),y)-V(F(f,2),y);
			auto sy = V(F(f,1),y)-V(F(f,2),y);
			double area = rx*sy - ry*sx;
			if(isnan(area))
				return 0.0;
			else
				return area;
		};

		dblA = DeriveddblA::Zero(m,1);
		for(size_t f = 0;f<m;f++)
		{
			for(int d = 0;d<3;d++)
			{
				const auto dblAd = proj_doublearea(d,(d+1)%3,f);
				//if(debugflag >12 && (dblAd > 0))
				//	cout<<dblAd<<"		";
				dblA(f) += dblAd*dblAd;

			}
		}
		dblA = dblA.array().sqrt().eval();debugflag++;
	}
}


MatrixXd V,U,tempU;
MatrixXi F,tempF;
SparseMatrix<double> L;
VectorXd area,interarea;
igl::opengl::glfw::Viewer viewer;
int sl = 2.0,vcount;
double areathreshold;
VectorXi EMAP;
MatrixXi E,EF,EI;

int main(int argc, char *argv[])
{
	// Load a mesh in OFF format
	igl::readOFF("mesh1.off", V, F);

	vcount = V.rows();
	U= V;

	/*const auto &key_down = [](igl::opengl::glfw::Viewer &viewer,unsigned char key,int mod)->bool
	{
		switch(key)
		{
			case ' ':
			{
				viewer.core.is_animating ^= 1;break;
			}
			default: return false;
		}
		return true;
	};

	const auto &pre_draw = [&](igl::opengl::glfw::Viewer & viewer) -> bool
	{
		if(viewer.core.is_animating)
		{*/
				//cout<<" ************ SPACE PRESSED *********** "<<endl;

				//L matrix
				igl::cotmatrix(V,F,L);

				//find average area for Wl
				igl::doublearea(V,F,area);
				area = area.array() / 2;
				double area_avg   = area.mean();

				cout<<"Total area of original mesh: "<<area.sum()<<endl;


				//RowVectorXd ringarea(vcount);
				ArrayXf ringarea = ArrayXf::Zero(vcount);

				//Define a one_ring area matrix for V
				for(int i = 0; i < F.rows(); i++)
				{
					ringarea(F(i,0)) += area(i);
					ringarea(F(i,1)) += area(i);
					ringarea(F(i,2)) += area(i);
				}


				//Create initial wl
				int i;
				MatrixXd wl(vcount,vcount);
				for(i=0;i<vcount;i++)
				{
					wl(i,i) = 0.001 * sqrt(area_avg);
				}

				//Initial wh
				MatrixXd wh(vcount,vcount);
				for(i=0;i<vcount;i++)
				{
					wh(i,i) = 1.0;
				}
				MatrixXd whi = wh;
				while(true)
				{
					MatrixXd a(vcount,vcount);
					a = wl * L;
					MatrixXd lhs(2*vcount,vcount);
					lhs<< a, wh;
					MatrixXd b(vcount,vcount);
					b = wh * U;
					ArrayXXd zro = ArrayXXd::Zero(vcount,3);
					MatrixXd rhs(2*vcount,3) ;
					rhs << zro, b;

					//Preserving the mesh incase the area threshold become 0
					tempU = U;
					tempF = F;

					//Solving for V(t+1)
					ColPivHouseholderQR<MatrixXd> solver(lhs);
					U= solver.solve(rhs);


				//viewer.data().clear();
				//viewer.data().set_mesh(U,F);
				//viewer.data().set_face_based(true);




					//Compare surface area of mesh with original
					mesh::doublearea(U,F,interarea);

					interarea = interarea.array() / 2;

					cout<<"Intermediate area: "<<interarea.sum()<<endl;

					areathreshold=interarea.sum()/area.sum();
					cout<<"Areathreshold: "<<areathreshold<<endl;
					cout<<"inter size: "<<interarea.size()<<endl;
					if(areathreshold < 0.001)
					{
						U = tempU;
						F = tempF;
						break;
					}
					igl::cotmatrix(U,F,L);
					//if(areathreshold < 0.3)
					//	cout<<L<<endl;

					//Update Wl for next iteration
					//Sl is suggested to be taken 2.0
					wl = wl * sl;

					//Calculate new one ring area
					ArrayXf interringarea = ArrayXf::Zero(vcount);
					for(int i = 0; i < F.rows(); i++)
					{
						interringarea(F(i,0)) += interarea(i);
						interringarea(F(i,1)) += interarea(i);
						interringarea(F(i,2)) += interarea(i);
					}

					//Update wh for next iteration
					for(i = 0;i < vcount; i++)
					{
						wh(i,i) = whi(i,i) * sqrt(ringarea(i)/interringarea(i));
					}

				}
				/*Second step:

				vector <MatrixXd> qmatrix;

				edge_flaps(F,E,EMAP,EF,EI);
				p = 0.5*(V.row(E(e,0))+V.row(E(e,1)));

				MatrixXd dummy = MatrixXd::Zero(4,4);
				for (int i = 1; i <= vcount; i++)
					qmatrix.push_back(dummy);

				for(i=0; i<F.rows(); i++)
				{
					for(int j=0; j<3; j++)
					{
						for(int add=1;add<=2;add++)
						{
							int k = (j+add)%3;
							MatrixXd kij = MatrixXd::Zero(3,4);
							VectorXd edgij(3);
							edgij << (U(F(i,k),0) - U(F(i,j),0)), (U(F(i,k),1) - U(F(i,j),1)), (U(F(i,j),2)-U(F(i,j),2));
							edgij = edgij.normalized();
							MatrixXd a = MatrixXd::Zero(3,3);
							a(0,1) = -edgij(2); a(1,0) = -a(0,1);
							a(0,2) = edgij(1); a(2,0) = -a(0,2);
							a(1,2) = -edgij(0); a(2,1) = -a(1,2);
							VectorXd b = -(edgij * U(F(i,j)));
							kij<<a,b;
							MatrixXd ktk = kij.transpose() * kij;
							qmatrix.at(F(i,j)) += ktk;
						}
					}
				}
				int dum=1;
				//When to stop: doubt
				while(dum<10)
				{
					MatrixXd normV = U.rowwise().homogeneous();
					VectorXd fa = VectorXd::Zero(vcount);
					for(i =0; i< vcount; i++)
					{
						VectorXd pt(4);
						pt << normV(i,0), normV(i,1), normV(i,2), normV(i,3);
						fa(i) = pt.transpose() * qmatrix.at(i) * pt;

					}
					VectorXd fbsum(vcount);
					for( i=0; i<F.rows(); i++)
					{
						for(int j=0; j<3; j++)
						{
							for(int add=1; add <=2; add++)
							{
								int k = (j+add)%3;
								double dist = sqrt( pow(( U(F(i,j),0) - U(F(i,k),0) ),2) +  pow(( U(F(i,j),1) - U(F(i,k),1) ),2) +  pow(( U(F(i,j),2) - U(F(i,k),2) ),2) );
								fbsum(F(i,j)) += dist;
							}
						}
					}
					double wa = 1.0, wb = 0.1;

					MatrixXd fmatrix(vcount,vcount);
					for(i = 0; i < vcount; i++)
					{
						for(int j=0; j< vcount; j++)
						{
							fmatrix(i,j) = -1;
						}
					}
					for(i =0; i < F.rows(); i++)
					{
						for(int j=0; j<3; j++)
						{
							for(int add=1; add<=2; add++)
							{
								int k=(j+add)%3;
								double dist =  sqrt( pow(( U(F(i,j),0) - U(F(i,k),0) ),2) +  pow(( U(F(i,j),1) - U(F(i,k),1) ),2) +  pow(( U(F(i,j),2) - U(F(i,k),2) ),2) );
								fmatrix(F(i,j),F(i,k)) = wa * ( fa(F(i,j)) + fa(F(i,k)) ) + wb * ( dist * fbsum(F(i,j)));
							}
						}
					}

					int meri=0,merj=0;
					for (i = 0; i< vcount; i++)
					{
						for (int j =0; j<vcount; j++)
						{
							if( ((fmatrix(i,j) > 0) && fmatrix(i,j) < fmatrix(meri,merj)) || fmatrix(meri, merj) < 0)
							{
								//cout<< i <<" "<<j<<"     "<<fmatrix(i,j)<<"   ################   "<<fmatrix(meri, merj)<<endl;
								meri = i; merj = j;
							}
						}
					}
					for(i = 0; i < F.rows(); i++)
					{
						for(int j = 0; j<3; j++)
						{
							if(F(i,j) == meri)
								F(i,j) = merj;
						}
					}
					U(meri) = U(merj);
					qmatrix.at(merj) += qmatrix.at(meri); dum++;
					cout<<"^^^^^^^^^^^^ Vertex "<<meri<<" merged to "<<merj<<" ^^^^^^^^^^^^^^"<<endl;
					viewer.data().set_vertices(U);
					viewer.core.align_camera_center(U,F);
				}*/
	//			break;
	//		}
	//
	//		default: return false;



		//}
		//return false;
	//};

	// Plot the mesh (Error: The new mesh has a different number of vertices/faces. Please clear the mesh before plotting.) doubt
	viewer.data().set_mesh(U, F);
	//viewer.callback_key_down = key_down;
	//viewer.data().set_face_based(true);
	//viewer.core.is_animating = true;
	//viewer.callback_key_down = key_down;
	//viewer.callback_pre_draw = pre_draw;
	return viewer.launch();

}
