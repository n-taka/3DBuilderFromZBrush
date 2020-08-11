#ifndef repairAs3DBuilder_CPP
#define repairAs3DBuilder_CPP

#include "repairAs3DBuilder.h"

#pragma warning(push)
#pragma warning(disable : 4018 4101 4129 4244 4267 4305 4566 4819 4996)
#include "igl/list_to_matrix.h"
#include "igl/polygon_mesh_to_triangle_mesh.h"
#include "igl/barycenter.h"
#include "igl/point_mesh_squared_distance.h"
#include "igl/barycentric_coordinates.h"
#include "igl/writeOBJ.h"
#pragma warning(pop)

#if defined _WIN64 | defined _WIN32
#pragma comment(lib, "windowsapp")
#include "winrt/Windows.Graphics.Printing3D.h"
#include "winrt/Windows.Storage.Streams.h"
#include "winrt/Windows.Storage.Pickers.h"
#endif

#include <iostream>

template <typename Scalar, typename Index>
void repairAs3DBuilder(
    const std::vector<std::vector<Scalar>> &V_in,
    const std::vector<std::vector<Index>> &F_in,
    std::vector<std::vector<Scalar>> &V_out,
    std::vector<std::vector<Index>> &F_out)
{
    const std::vector<std::vector<std::pair<Scalar, Scalar>>> UV_in;
    const std::vector<std::vector<Scalar>> VC_in;
    const std::vector<Scalar> M_in;
    const std::vector<Index> G_in;
    std::vector<std::vector<std::pair<Scalar, Scalar>>> UV_out;
    std::vector<std::vector<Scalar>> VC_out;
    std::vector<Scalar> M_out;
    std::vector<Index> G_out;

    repairAs3DBuilder(V_in, F_in, UV_in, VC_in, M_in, G_in,
        V_out, F_out, UV_out, VC_out, M_out, G_out);
}

template <typename Scalar, typename Index>
void repairAs3DBuilder(
    const std::vector<std::vector<Scalar>> &V_in,
    const std::vector<std::vector<Index>> &F_in,
    const std::vector<std::vector<std::pair<Scalar, Scalar>>> &UV_in,
    const std::vector<std::vector<Scalar>> &VC_in,
    const std::vector<Scalar> &M_in,
    const std::vector<Index> &G_in,
    std::vector<std::vector<Scalar>> &V_out,
    std::vector<std::vector<Index>> &F_out,
    std::vector<std::vector<std::pair<Scalar, Scalar>>> &UV_out,
    std::vector<std::vector<Scalar>> &VC_out,
    std::vector<Scalar> &M_out,
    std::vector<Index> &G_out)
{
    #if defined _WIN64 | defined _WIN32
    // convert to 3MF
    winrt::Windows::Graphics::Printing3D::Printing3DModel model3MF;
    model3MF.Unit(winrt::Windows::Graphics::Printing3D::Printing3DModelUnit::Millimeter);
    model3MF.Metadata().Insert(L"Title", L"Repair");
    model3MF.Metadata().Insert(L"Designer", L"Kazutaka Nakashima");

    // prepare mesh
    winrt::Windows::Graphics::Printing3D::Printing3DMesh mesh3MF;
    // vertices
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> V;
    {
        Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> Vs;
        igl::list_to_matrix(V_in, Vs);
        V = Vs.template cast<double>();

        winrt::Windows::Graphics::Printing3D::Printing3DBufferDescription description;
        description.Format = winrt::Windows::Graphics::Printing3D::Printing3DBufferFormat::Printing3DDouble;
        mesh3MF.VertexPositionsDescription(description);
        // have 3 xyz values
        description.Stride = 3;
        // update vertex count
        mesh3MF.VertexCount(V.rows());
        // write data
        mesh3MF.CreateVertexPositions(sizeof(double) * V.size());
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> VRowMajor = V;

        winrt::Windows::Storage::Streams::IBuffer buffer = mesh3MF.GetVertexPositions();
        std::memcpy(buffer.data(), VRowMajor.data(), VRowMajor.size() * sizeof(double));
    }
    // faces
    // 3MF format assumes pure triangle mesh.
    Eigen::Matrix<Index, Eigen::Dynamic, Eigen::Dynamic> F;
    std::vector<std::vector<std::pair<Scalar, Scalar>>> UV;
    std::vector<Index> G;
    {
        // triangulate
        {
            int faceCount = 0;
            for (const auto& f : F_in)
            {
                faceCount += (f.size() -2);
            }
            F.resize(faceCount, 3);
            if (UV_in.size() > 0)
            {
                UV.resize(faceCount, std::vector<std::pair<Scalar, Scalar>>(3));
            }
            if (G_in.size() > 0)
            {
                G.resize(faceCount);
            }
            int faceIdx = 0;
            for (int fIdx=0;fIdx<F_in.size();++fIdx)
            {
                const std::vector<Index> &f = F_in.at(fIdx);
                for (int tri=0;tri<f.size()-2;++tri)
                {
                    F.row(faceIdx) << f.at(0), f.at(tri+1), f.at(tri+2);
                    if (UV_in.size() > 0)
                    {
                        UV.at(faceIdx).at(0) = UV_in.at(fIdx).at(0);
                        UV.at(faceIdx).at(1) = UV_in.at(fIdx).at(tri+1);
                        UV.at(faceIdx).at(2) = UV_in.at(fIdx).at(tri+2);
                    }
                    if (G_in.size() > 0)
                    {
                        G.at(faceIdx) = G_in.at(fIdx);
                    }
                    faceIdx++;
                }
            }
        }

        winrt::Windows::Graphics::Printing3D::Printing3DBufferDescription description;
        description.Format = winrt::Windows::Graphics::Printing3D::Printing3DBufferFormat::Printing3DUInt;
        mesh3MF.TriangleIndicesDescription(description);
        // triangle
        description.Stride = 3;
        mesh3MF.IndexCount(F.rows());

        mesh3MF.CreateTriangleIndices(sizeof(unsigned int) * F.size());
        Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> FRowMajor = F.template cast<unsigned int>();
        ;
        winrt::Windows::Storage::Streams::IBuffer buffer = mesh3MF.GetTriangleIndices();

        std::memcpy(buffer.data(), FRowMajor.data(), FRowMajor.size() * sizeof(unsigned int));
    }

    // append to model
    model3MF.Meshes().Append(mesh3MF);

    // perform repair
    winrt::Windows::Foundation::IAsyncAction repair = model3MF.RepairAsync();
    repair.get();

    winrt::Windows::Graphics::Printing3D::Printing3DMesh repairedMesh3MF;
    repairedMesh3MF = model3MF.Meshes().GetAt(0);

    // For projecting attributes such as mask and vertex color,
    // we first get in Eigen Matrix style and then convert to std::vector
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> V_repair;
    Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> F_repair;
    // vertices
    {
        winrt::Windows::Storage::Streams::IBuffer buffer = repairedMesh3MF.GetVertexPositions();

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> VRowMajor;
        VRowMajor.resize(repairedMesh3MF.VertexCount(), 3);

        std::memcpy(VRowMajor.data(), buffer.data(), VRowMajor.size() * sizeof(double));

        V_repair = VRowMajor;
    }
    // faces
    {
        winrt::Windows::Storage::Streams::IBuffer buffer = repairedMesh3MF.GetTriangleIndices();

        Eigen::Matrix<unsigned int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> FRowMajor;
        FRowMajor.resize(repairedMesh3MF.IndexCount(), 3);

        std::memcpy(FRowMajor.data(), buffer.data(), FRowMajor.size() * sizeof(unsigned int));

        F_repair = FRowMajor.template cast<int>();
    }

    // write the result
    // Vertices
    V_out.resize(V_repair.rows(), std::vector<Scalar>(V_repair.cols()));
    for (int v=0;v<V_repair.rows();++v)
    {
        for (int xyz=0;xyz<V_repair.cols();++xyz)
        {
            V_out.at(v).at(xyz) = static_cast<Scalar>(V_repair(v, xyz));
        }
    }
    // Faces
    F_out.resize(F_repair.rows(), std::vector<Index>(F_repair.cols()));
    for (int f=0;f<F_repair.rows();++f)
    {
        for (int xyz=0;xyz<F_repair.cols();++xyz)
        {
            F_out.at(f).at(xyz) = static_cast<Index>(F_repair(f, xyz));
        }
    }

    // calculate nearest point for mapping attributes
    const bool needVertAttrib = (UV_in.size() > 0 || VC_in.size() > 0 || M_in.size() > 0);
    const bool needFaceAttrib = (G_in.size() > 0);
    Eigen::Matrix<int, Eigen::Dynamic, 1> vertI, faceI;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vertBC, faceBC;
    if (needVertAttrib || needFaceAttrib)
    {
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> P;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vertP;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> faceP;
        if (needVertAttrib)
        {
            vertP = V_repair;
        }
        else
        {
            vertP.resize(0, 3);
        }
        if (needFaceAttrib)
        {
            igl::barycenter(V_repair, F_repair, faceP);
        }
        else
        {
            faceP.resize(0, 3);
        }
        P.resize(vertP.rows() + faceP.rows(), 3);
        P << vertP,
            faceP;

        Eigen::Matrix<double, Eigen::Dynamic, 1> sqrD;
        Eigen::Matrix<int, Eigen::Dynamic, 1> I;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> C;
        igl::point_mesh_squared_distance(P, V, F, sqrD, I, C);

        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> vertC, faceC;
        if (needVertAttrib)
        {
            vertI = I.block(0, 0, vertP.rows(), 1);
            vertC = C.block(0, 0, vertP.rows(), 3);
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A_, B_, C_;
            A_.resize(vertI.rows(), 3);
            B_.resize(vertI.rows(), 3);
            C_.resize(vertI.rows(), 3);
            for (int v=0;v<vertI.rows();++v)
            {
                A_.row(v) = V.row(F(vertI(v, 0), 0));
                B_.row(v) = V.row(F(vertI(v, 0), 1));
                C_.row(v) = V.row(F(vertI(v, 0), 2));
            }
            igl::barycentric_coordinates(vertC, A_, B_, C_, vertBC);
        }
        else
        {
            vertI.resize(0, 1);
            vertC.resize(0, 3);
            vertBC.resize(0, 3);
        }
        if (needFaceAttrib)
        {
            faceI = I.block(vertP.rows(), 0, faceP.rows(), 1);
            faceC = C.block(vertP.rows(), 0, faceP.rows(), 3);
            Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> A_, B_, C_;
            A_.resize(faceI.rows(), 3);
            B_.resize(faceI.rows(), 3);
            C_.resize(faceI.rows(), 3);
            for (int v=0;v<faceI.rows();++v)
            {
                A_.row(v) = V.row(F(faceI(v, 0), 0));
                B_.row(v) = V.row(F(faceI(v, 0), 1));
                C_.row(v) = V.row(F(faceI(v, 0), 2));
            }
            igl::barycentric_coordinates(faceC, A_, B_, C_, faceBC);
        }
        else
        {
            faceI.resize(0, 1);
            faceC.resize(0, 3);
            faceBC.resize(0, 3);
        }
    }


    UV_out.clear();
    VC_out.clear();
    M_out.clear();
    G_out.clear();
    // UV
    if (UV_in.size() > 0)
    {
        UV_out.resize(F_repair.rows(), std::vector<std::pair<Scalar, Scalar>>(F_repair.cols()));
    }

    // VC
    if (VC_in.size() > 0)
    {
        VC_out.resize(V_repair.rows(), std::vector<Scalar>(4));
    }

    // M
    if (M_in.size() > 0)
    {
        M_out.resize(V_repair.rows());
    }

    // G
    if (G_in.size() > 0)
    {
        G_out.resize(F_repair.rows());
    }

    for (int f=0;f<F_repair.rows();++f)
    {
        for (int fv=0;fv<F_repair.cols();++fv)
        {
            // UV
            if (UV_in.size() > 0)
            {
                UV_out.at(f).at(fv).first = static_cast<Scalar>(0.0);
                UV_out.at(f).at(fv).second = static_cast<Scalar>(0.0);
                const int vIdx = F_repair(f, fv);
                const int fIdxInOrig = vertI(vIdx, 0);
                for (int uvw=0;uvw<vertBC.cols();++uvw)
                {
                    UV_out.at(f).at(fv).first += (vertBC(vIdx, uvw) * UV.at(fIdxInOrig).at(uvw).first);
                    UV_out.at(f).at(fv).second += (vertBC(vIdx, uvw) * UV.at(fIdxInOrig).at(uvw).second);
                }
            }
            // G
            if (G_in.size() > 0)
            {
                const int fIdxInOrig = faceI(f, 0);
                G_out.at(f) = G.at(fIdxInOrig);
            }
        }
    }

    for (int v=0;v<V_repair.rows();++v)
    {
        // M
        if (M_in.size() > 0)
        {
            M_out.at(v) = static_cast<Scalar>(0.0);
            const int fIdxInOrig = vertI(v, 0);
            for (int uvw=0;uvw<vertBC.cols();++uvw)
            {
                M_out.at(v) += (static_cast<Scalar>(vertBC(v, uvw)) * M_in.at(F(fIdxInOrig, uvw)));
            }
        }

        // VC
        if (VC_in.size() > 0)
        {
            for (int rgba=0;rgba<4;++rgba)
            {
                VC_out.at(v).at(rgba) = static_cast<Scalar>(0.0);
            }
            const int fIdxInOrig = vertI(v, 0);
            for (int uvw=0;uvw<vertBC.cols();++uvw)
            {
                for (int rgba=0;rgba<4;++rgba)
                {
                    VC_out.at(v).at(rgba) += (vertBC(v, uvw) * VC_in.at(F(fIdxInOrig, uvw)).at(rgba));
                }
            }
        }
    }

    #endif
}

#endif
