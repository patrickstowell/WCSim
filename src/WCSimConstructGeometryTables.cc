#include "WCSimDetectorConstruction.hh"
#include "WCSimPmtInfo.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4Vector3D.hh"
#include "globals.hh"
#include "G4VisAttributes.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4VPVParameterisation.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include <sstream>
#include <iomanip>

#include <fstream>
#include <vector>
#include <string>


using std::setw;
// These routines are object registration routines that you can pass
// to the traversal code.

void WCSimDetectorConstruction::PrintGeometryTree
(G4VPhysicalVolume* aPV ,int aDepth, int /*replicaNo*/, 
 const G4Transform3D& aTransform) 
{
  for (int levels = 0; levels < aDepth; levels++) G4cout << " ";
  G4cout << aPV->GetName() << " Level:" << aDepth 
	 << " Pos:" << aTransform.getTranslation() 
	 << " Rot:" << aTransform.getRotation().getTheta()/deg 
	 << "," << aTransform.getRotation().getPhi()/deg 
	 << "," << aTransform.getRotation().getPsi()/deg
	 << G4endl;
}

void WCSimDetectorConstruction::GetWCGeom
(G4VPhysicalVolume* aPV ,int aDepth, int /*replicaNo*/, 
 const G4Transform3D& aTransform) 
{

    // Grab mishmash of useful information from the tree while traversing it
    // This information will later be written to the geometry file
    // (Alternatively one might define accessible constants)
  
    if ((aPV->GetName() == "WCBox")) {    // last condition is the HyperK Envelope name.
    // Stash info in data member
        WCOffset = G4ThreeVector(aTransform.getTranslation().getX()/cm,
			     aTransform.getTranslation().getY()/cm,
			     aTransform.getTranslation().getZ()/cm);
        WCXRotation = G4ThreeVector(aTransform.getRotation().xx(), aTransform.getRotation().xy(), aTransform.getRotation().xz());
        WCYRotation = G4ThreeVector(aTransform.getRotation().yx(), aTransform.getRotation().yy(), aTransform.getRotation().yz());
        WCZRotation = G4ThreeVector(aTransform.getRotation().zx(), aTransform.getRotation().zy(), aTransform.getRotation().zz());
    }

    // Stash info in data member
    // AH Need to store this in CM for it to be understood by SK code
    WCPMTSize = WCPMTRadius/cm;// I think this is just a variable no if needed
    WCPMTSize2 = WCPMTRadius2/cm;// I think this is just a variable no if needed

    // Note WC can be off-center... get both extremities
    static G4double zmin=100000,zmax=-100000.;
    static G4double xmin=100000,xmax=-100000.;
    static G4double ymin=100000,ymax=-100000.;
    if (aDepth == 0) { // Reset for this traversal
        xmin=100000,xmax=-100000.; 
        ymin=100000,ymax=-100000.; 
        zmin=100000,zmax=-100000.; 
    }

    if ((aPV->GetName() == "WCCapBlackSheet") || (aPV->GetName().find("glassFaceWCPMT") != std::string::npos)){ 
      G4double x =  aTransform.getTranslation().getX()/cm;
      G4double y =  aTransform.getTranslation().getY()/cm;
      G4double z =  aTransform.getTranslation().getZ()/cm;
      
      if (x<xmin){xmin=x;}
      if (x>xmax){xmax=x;}
      
      if (y<ymin){ymin=y;}
      if (y>ymax){ymax=y;}

      if (z<zmin){zmin=z;}
      if (z>zmax){zmax=z;}
      
      WCCylInfo[0] = xmax-xmin;
      WCCylInfo[1] = ymax-ymin;
      WCCylInfo[2] = zmax-zmin;
      //      G4cout << "determine height: " << zmin << "  " << zmax << " " << aPV->GetName()<<" " << z  << G4endl;
  } 
}



void WCSimDetectorConstruction::DescribeAndRegisterPMT(G4VPhysicalVolume* aPV ,int aDepth, int replicaNo,
                                                       const G4Transform3D& aTransform) 
{
  static std::string replicaNoString[20];

  std::stringstream depth;
  std::stringstream pvname;

  depth << replicaNo;
  pvname << aPV->GetName();

  replicaNoString[aDepth] = pvname.str() + "-" + depth.str();

 
  //TF: To Consider: add a separate table for mPMT positions? Need to use its orientation anyway
  // Could be useful for the near future. Need to add an == WCMultiPMT here then.
  if (aPV->GetName()== WCIDCollectionName || aPV->GetName()== WCIDCollectionName2 ||aPV->GetName()== WCODCollectionName ) 
    {

    // First increment the number of PMTs in the tank.
      if(aPV->GetName()== WCIDCollectionName) totalNumPMTs++;
      else if(aPV->GetName()== WCIDCollectionName2) totalNumPMTs2++;
      else if(aPV->GetName()== WCODCollectionName) totalNumODPMTs++;

    // Put the location of this tube into the location map so we can find
    // its ID later.  It is coded by its tubeTag string.
    // This scheme must match that used in WCSimWCSD::ProcessHits()

    std::string tubeTag;
    int mPMT_pmtno = -1;
    bool foundString = false;
    int mPMTid = -1;
    //
    
    for (int i=0; i <= aDepth; i++){
      tubeTag += ":" + replicaNoString[i];
      std::string::size_type position = replicaNoString[i].find("pmt-");
      if( position == 0) {
        foundString = true;
        mPMT_pmtno = atoi(replicaNoString[i].substr(position+4).c_str())+1;
        if(mPMT_pmtno == 1) {
          if(aPV->GetName()== WCIDCollectionName) totalNum_mPMTs++;
          else if(aPV->GetName()== WCIDCollectionName2) totalNum_mPMTs2++;
        }
      }
      position = replicaNoString[i].find("WCMultiPMT-");
      if(position == 0) mPMTid = atoi(replicaNoString[i].substr(position+11).c_str());
    }
    if(!foundString){
      // to distinguish mPMT PMTs from single PMTs:
      mPMT_pmtno = 0;
      if(aPV->GetName()== WCIDCollectionName) totalNum_mPMTs++;
      else if(aPV->GetName()== WCIDCollectionName2) totalNum_mPMTs2++;
    }
    
    // G4cout << tubeTag << G4endl;
    if(aPV->GetName()== WCIDCollectionName){    
      if ( tubeLocationMap.find(tubeTag) != tubeLocationMap.end() ) {
        G4cerr << "Repeated tube tag: " << tubeTag << G4endl;
        G4cerr << "Assigned to both tube #" << tubeLocationMap[tubeTag] << " and #" << totalNumPMTs << G4endl;
        G4cerr << "Cannot continue -- hits will not be recorded correctly." << G4endl;
        G4cerr << "Please make sure that logical volumes with multiple placements are each given a unique copy number"
               << G4endl;
        assert(false);
      }
      tubeLocationMap[tubeTag] = totalNumPMTs;
      // Put the transform for this tube into the map keyed by its ID
      tubeIDMap[totalNumPMTs] = aTransform;

      if (!useReplica&&readFromTable) mPMTIDMap[totalNumPMTs] = std::make_pair(mPMTid,mPMT_pmtno);
      else mPMTIDMap[totalNumPMTs] = std::make_pair(totalNum_mPMTs,mPMT_pmtno);
    }//ID PMT 1

    else if(aPV->GetName()== WCIDCollectionName2){    
      if ( tubeLocationMap2.find(tubeTag) != tubeLocationMap2.end() ) {
        G4cerr << "Repeated tube tag: " << tubeTag << G4endl;
        G4cerr << "Assigned to both tube #" << tubeLocationMap2[tubeTag] << " and #" << totalNumPMTs2 << G4endl;
        G4cerr << "Cannot continue -- hits will not be recorded correctly."  << G4endl;
        G4cerr << "Please make sure that logical volumes with multiple placements are each given a unique copy number" << G4endl;
        assert(false);
      }
      tubeLocationMap2[tubeTag] = totalNumPMTs2;
      // Put the transform for this tube into the map keyed by its ID
      tubeIDMap2[totalNumPMTs2] = aTransform;

      if (!useReplica&&readFromTable) mPMTIDMap2[totalNumPMTs2] = std::make_pair(mPMTid,mPMT_pmtno);
      mPMTIDMap2[totalNumPMTs2] = std::make_pair(totalNum_mPMTs2,mPMT_pmtno);
    }//ID PMT 2

    else if(aPV->GetName()== WCODCollectionName) {
      if(ODtubeLocationMap.find(tubeTag) != ODtubeLocationMap.end()) {
        G4cerr << "Repeated tube tag: " << tubeTag << G4endl;
        G4cerr << "Assigned to both tube #" << ODtubeLocationMap[tubeTag] << " and #" << totalNumODPMTs << G4endl;
        G4cerr << "Cannot continue -- hits will not be recorded correctly." << G4endl;
        G4cerr << "Please make sure that logical volumes with multiple placements are each given a unique copy number"
               << G4endl;
        assert(false);
      }
      ODtubeLocationMap[tubeTag] = totalNumODPMTs;

      // Put the transform for this tube into the map keyed by its ID
      ODtubeIDMap[totalNumODPMTs] = aTransform;

      mPMTODMap[totalNumODPMTs] = std::make_pair(totalNumODPMTs,mPMT_pmtno);
    }//OD PMT

    // G4cout <<  "depth " << depth.str() << G4endl;
    // G4cout << "tubeLocationmap[" << tubeTag  << "]= " << tubeLocationMap[tubeTag] << "\n";

    // Print
    //     G4cout << "Tube: "<<std::setw(4) << totalNumPMTs << " " << tubeTag
    //	   << " Pos:" << aTransform.getTranslation()/cm 
    //	   << " Rot:" << aTransform.getRotation().getTheta()/deg 
    //	   << "," << aTransform.getRotation().getPhi()/deg 
    //	   << "," << aTransform.getRotation().getPsi()/deg
    //	   << G4endl; 
    }//if ID PMT 1/2 or OD PMT
}

// Utilities to do stuff with the info we have found.


std::vector< std::vector<int> > ObtainSTLFacets(const G4Polyhedron& obj){
    std::vector<std::vector<int>> data;

    G4int iFace;
    G4int n; 
    G4int iNodes[100];
    G4int edgeFlags[100];
    G4int iFaces[100];
    
    for (iFace = 0; iFace < obj.GetNoFacets(); iFace++) {
        obj.GetFacet(iFace+1, n, iNodes, edgeFlags, iFaces);
        std::vector<int> temp;
        if (n == 4){
            temp.push_back(iFaces[0] - 1);
            temp.push_back(edgeFlags[0]-1);
            temp.push_back(iNodes[0] - 1);
            temp.push_back(iNodes[1] - 1);
            temp.push_back(iNodes[2] - 1);
            temp.push_back(iNodes[3] - 1);        
        } else {
            temp.push_back(iFaces[0] - 1);
            temp.push_back(edgeFlags[0]-1);
            temp.push_back(iNodes[0] - 1);
            temp.push_back(iNodes[1] - 1);
            temp.push_back(iNodes[2] - 1);
            temp.push_back(iNodes[0] - 1); 
        }
        data.push_back(temp);
    }    
    return data;
}


void ExportToPLY(const std::string& filename,
                 const std::vector<std::vector<double>>& trans_vertices,
                 const std::vector<std::vector<int>>& trans_facets)
{
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "Could not open file for writing: " << filename << std::endl;
        return;
    }

    size_t num_vertices = trans_vertices.size();
    size_t num_faces = trans_facets.size();

    out << "ply\n";
    out << "format ascii 1.0\n";
    out << "element vertex " << num_vertices << "\n";
    out << "property float x\n";
    out << "property float y\n";
    out << "property float z\n";
    out << "element face " << num_faces << "\n";
    out << "property list uchar int vertex_index\n";
    out << "end_header\n";

    // Write vertices
    for (const auto& v : trans_vertices) {
        out << v[0] << " " << v[1] << " " << v[2] << "\n";
    }

    // Write faces (triangles)
    for (const auto& f : trans_facets) {
        out << "3 " << static_cast<int>(f[0]) << " " << static_cast<int>(f[1]) << " " << static_cast<int>(f[2]) << "\n";
    }

    out.close();
    std::cout << "PLY export complete: " << filename << std::endl;
}





#include <fstream>
#include <vector>
#include <string>
#include <iostream>

void ExportToOBJ_Append(const std::string& filename,
                        const std::vector<std::vector<double>>& vertices,
                        const std::vector<std::vector<int>>& faces,
                        const std::string& object_name,
                        size_t& global_vertex_offset)
{
    std::ofstream out(filename, std::ios::app);  // append mode
    if (!out.is_open()) {
        std::cerr << "Failed to open OBJ file for appending: " << filename << std::endl;
        return;
    }

    // Write object name
    out << "o " << object_name << "\n";

    // Write vertices
    for (const auto& v : vertices) {
        out << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
    }

    // Write faces with correct global offset
    for (const auto& f : faces) {
        out << "f "
            << (f[0] + global_vertex_offset) << " "
            << (f[1] + global_vertex_offset) << " "
            << (f[2] + global_vertex_offset) << "\n";
    }

    // Update the global offset for the next object
    global_vertex_offset += vertices.size();

    out.close();
}


// void ExportToOBJ(const std::string& filename,
//                  const std::vector<std::vector<std::vector<double>>>& all_vertices,
//                  const std::vector<std::vector<std::vector<double>>>& all_faces,
//                  const std::vector<std::string>& object_names = {})
// {
//     std::ofstream out(filename);
//     if (!out.is_open()) {
//         std::cerr << "Failed to open OBJ file for writing: " << filename << std::endl;
//         return;
//     }

//     size_t global_vertex_offset = 1;  // OBJ indices are 1-based

//     for (size_t obj_idx = 0; obj_idx < all_vertices.size(); ++obj_idx) {
//         const auto& vertices = all_vertices[obj_idx];
//         const auto& faces = all_faces[obj_idx];
//         std::string obj_name = (object_names.size() > obj_idx) ? object_names[obj_idx] : "Object" + std::to_string(obj_idx);

//         // Write object name
//         out << "o " << obj_name << "\n";

//         // Write vertices
//         for (const auto& v : vertices) {
//             out << "v " << v[0] << " " << v[1] << " " << v[2] << "\n";
//         }

//         // Write faces (add offset to local indices)
//         for (const auto& f : faces) {
//             out << "f " 
//                 << (f[0] + global_vertex_offset) << " "
//                 << (f[1] + global_vertex_offset) << " "
//                 << (f[2] + global_vertex_offset) << "\n";
//         }

//         global_vertex_offset += vertices.size();
//     }

//     out.close();
//     std::cout << "OBJ export complete: " << filename << std::endl;
// }


void WCSimDetectorConstruction::BuildSTLModel(G4VPhysicalVolume* aPV ,int aDepth, int replicaNo,
                                                       const G4Transform3D& aTransform) 
{
  if (aDepth != 5) return;
  // if (global_obj_size > 40000) return;

  static std::string replicaNoString[20];

  std::stringstream depth;
  std::stringstream pvname;
  std::stringsstream filename_in;
  filename_in << "plylayer_" << aDepth;
  

  depth << replicaNo << "-" << aDepth;
  pvname << aPV->GetName();
  
  if (std::find(global_obj_size.begin(), global_obj_size.end(), aDepth) == global_obj_size.end()){
    global_obj_size[aDepth] = 1;
  }
  

 // Check material outer is water, if so skip it.
  replicaNoString[aDepth] = pvname.str() + "-" + depth.str();

  auto poly = aPV->GetLogicalVolume()->GetSolid()->CreatePolyhedron();

  std::vector< std::vector<double> > vertices;
  for (int i = 0; i < poly->GetNoVertices(); i++){
    auto v = poly->GetVertex(i+1);
    vertices.push_back( {v[0], v[1], v[2]} );
  }
  
  auto facets = ObtainSTLFacets(*poly);

  // std::vector< std::vector<double> > normals;
  // for (int i = 0; i < poly->GetNoFacets(); i++){
  //   auto v = poly->GetUnitNormal(i+1);
  //   normals.push_back( {v[0], v[1], v[2]} );
  // }
  
  std::vector< std::vector<double> > trans_vertices;
  for (auto const v : vertices){
    auto p = G4ThreeVector(v[0], v[1], v[2]);
    p = aTransform.getRotation()*p + aTransform.getTranslation();
    trans_vertices.push_back( {p[0], p[1], p[2]} );
  }

  // std::vector< std::vector<double> > trans_normals;
  // for (auto const n : normals){
  //   auto p = G4ThreeVector(n[0], n[1], n[2]);
  //   p = aTransform.getRotation()*p;
  //   trans_normals.push_back( {p[0], p[1], p[2]} );
  // }

  std::vector< std::vector<int> > trans_facets;
  for (auto const f : facets){
    std::vector<int> ff = {f[2], f[3], f[4], f[5]};

    trans_facets.push_back(std::vector<int>({ff[0], ff[1], ff[2]}));
    trans_facets.push_back(std::vector<int>({ff[0], ff[1], ff[3]}));
    trans_facets.push_back(std::vector<int>({ff[0], ff[2], ff[1]}));
    trans_facets.push_back(std::vector<int>({ff[0], ff[2], ff[3]}));
    trans_facets.push_back(std::vector<int>({ff[0], ff[3], ff[2]}));
    trans_facets.push_back(std::vector<int>({ff[0], ff[3], ff[3]}));

    trans_facets.push_back(std::vector<int>({ff[1], ff[0], ff[2]}));
    trans_facets.push_back(std::vector<int>({ff[1], ff[0], ff[3]}));
    trans_facets.push_back(std::vector<int>({ff[1], ff[2], ff[0]}));
    trans_facets.push_back(std::vector<int>({ff[1], ff[2], ff[3]}));
    trans_facets.push_back(std::vector<int>({ff[1], ff[3], ff[0]}));
    trans_facets.push_back(std::vector<int>({ff[1], ff[3], ff[3]}));

    trans_facets.push_back(std::vector<int>({ff[2], ff[0], ff[1]}));
    trans_facets.push_back(std::vector<int>({ff[2], ff[0], ff[3]}));
    trans_facets.push_back(std::vector<int>({ff[2], ff[1], ff[0]}));
    trans_facets.push_back(std::vector<int>({ff[2], ff[1], ff[3]}));
    trans_facets.push_back(std::vector<int>({ff[2], ff[3], ff[0]}));
    trans_facets.push_back(std::vector<int>({ff[2], ff[3], ff[1]}));
  }

  std::map<int, int> local_global_vertex_map;

  std::vector< std::vector<double> > group_vertices;
  for (int i = 0; i < trans_vertices.size(); i++){
    int local_id = i;
    int global_id = i; //global_vertices.size();
    local_global_vertex_map[local_id] = global_id;
    // global_vertices.push_back( trans_vertices[i] );
    group_vertices.push_back( trans_vertices[i] );
  }    

  std::vector< std::vector<int> > group_facets;

  for (int i = 0; i < trans_facets.size(); i++){

    auto f_local = trans_facets[i];
    std::vector<int> f_global;
    f_global.push_back(local_global_vertex_map[f_local[0]]);
    f_global.push_back(local_global_vertex_map[f_local[1]]);
    f_global.push_back(local_global_vertex_map[f_local[2]]);
    // global_facets.push_back( f_global );
    group_facets.push_back( f_global );

  }

  std::cout << "FACET SIZE : " << aDepth << " " << global_obj_size[aDepth] << std::endl;

  ExportToOBJ_Append(filename_in.str(), group_vertices, group_facets, replicaNoString[aDepth], global_obj_size[aDepth]);     

}

void WCSimDetectorConstruction::DumpSTLToFile(){
  // ExportToPLY("test.ply", global_vertices, global_facets);
  
}




// void WCSimDetectorConstruction::BuildSTLModel(G4VPhysicalVolume* aPV ,int aDepth, int replicaNo,
//                                                        const G4Transform3D& aTransform) 
// {
//   if (replicaNo > 1) return;
//   if (global_facets.size() > 400000000) return;

//   static std::string replicaNoString[20];

//   std::stringstream depth;
//   std::stringstream pvname;

//   depth << replicaNo;
//   pvname << aPV->GetName();


//   replicaNoString[aDepth] = pvname.str() + "-" + depth.str();

//   auto poly = aPV->GetLogicalVolume()->GetSolid()->CreatePolyhedron();

//   std::vector< std::vector<double> > vertices;
//   for (int i = 0; i < poly->GetNoVertices(); i++){
//     auto v = poly->GetVertex(i+1);
//     vertices.push_back( {v[0], v[1], v[2]} );
//   }
  
//   auto facets = ObtainSTLFacets(*poly);

//   // std::vector< std::vector<double> > normals;
//   // for (int i = 0; i < poly->GetNoFacets(); i++){
//   //   auto v = poly->GetUnitNormal(i+1);
//   //   normals.push_back( {v[0], v[1], v[2]} );
//   // }
  
//   std::vector< std::vector<double> > trans_vertices;
//   for (auto const v : vertices){
//     auto p = G4ThreeVector(v[0], v[1], v[2]);
//     p = aTransform.getRotation()*p + aTransform.getTranslation();
//     trans_vertices.push_back( {p[0], p[1], p[2]} );
//   }

//   // std::vector< std::vector<double> > trans_normals;
//   // for (auto const n : normals){
//   //   auto p = G4ThreeVector(n[0], n[1], n[2]);
//   //   p = aTransform.getRotation()*p;
//   //   trans_normals.push_back( {p[0], p[1], p[2]} );
//   // }

//   std::vector< std::vector<int> > trans_facets;
//   for (auto const f : facets){
//     std::vector<int> ff = {f[2], f[3], f[4], f[5]};

//     trans_facets.push_back(std::vector<int>({ff[0], ff[1], ff[2]}));
//     trans_facets.push_back(std::vector<int>({ff[0], ff[1], ff[3]}));
//     trans_facets.push_back(std::vector<int>({ff[0], ff[2], ff[1]}));
//     trans_facets.push_back(std::vector<int>({ff[0], ff[2], ff[3]}));
//     trans_facets.push_back(std::vector<int>({ff[0], ff[3], ff[2]}));
//     trans_facets.push_back(std::vector<int>({ff[0], ff[3], ff[3]}));

//     trans_facets.push_back(std::vector<int>({ff[1], ff[0], ff[2]}));
//     trans_facets.push_back(std::vector<int>({ff[1], ff[0], ff[3]}));
//     trans_facets.push_back(std::vector<int>({ff[1], ff[2], ff[0]}));
//     trans_facets.push_back(std::vector<int>({ff[1], ff[2], ff[3]}));
//     trans_facets.push_back(std::vector<int>({ff[1], ff[3], ff[0]}));
//     trans_facets.push_back(std::vector<int>({ff[1], ff[3], ff[3]}));

//     trans_facets.push_back(std::vector<int>({ff[2], ff[0], ff[1]}));
//     trans_facets.push_back(std::vector<int>({ff[2], ff[0], ff[3]}));
//     trans_facets.push_back(std::vector<int>({ff[2], ff[1], ff[0]}));
//     trans_facets.push_back(std::vector<int>({ff[2], ff[1], ff[3]}));
//     trans_facets.push_back(std::vector<int>({ff[2], ff[3], ff[0]}));
//     trans_facets.push_back(std::vector<int>({ff[2], ff[3], ff[1]}));
//   }

//   std::map<int, int> local_global_vertex_map;
//   for (int i = 0; i < trans_vertices.size(); i++){
//     int local_id = i;
//     int global_id = global_vertices.size();
//     local_global_vertex_map[local_id] = global_id;
//     global_vertices.push_back( trans_vertices[i] );
//   }    

//   for (int i = 0; i < trans_facets.size(); i++){

//     auto f_local = trans_facets[i];
//     std::vector<int> f_global;
//     f_global.push_back(local_global_vertex_map[f_local[0]]);
//     f_global.push_back(local_global_vertex_map[f_local[1]]);
//     f_global.push_back(local_global_vertex_map[f_local[2]]);
//     global_facets.push_back( f_global );

//   }

//   std::cout << "FACET SIZE : " << global_facets.size() << std::endl;

                
  
// }

// void WCSimDetectorConstruction::DumpSTLToFile(){
//   ExportToPLY("test.ply", global_vertices, global_facets);
  
// }


// Output to WC geometry text file
void WCSimDetectorConstruction::DumpGeometryTableToFile()
{
  // Open a file
  std::string filename = "geofile_" + WCDetectorName + ".txt";
  //geoFile.open("geofile.txt", std::ios::out);
  geoFile.open(filename.c_str(), std::ios::out);

  geoFile.precision(2);
  geoFile.setf(std::ios::fixed);

  // (JF) Get first tube transform for filling in detector radius
  // the height is still done with WCCylInfo above
  G4Transform3D firstTransform = tubeIDMap[2];
  innerradius = sqrt(pow(firstTransform.getTranslation().getX()/cm,2)
                            + pow(firstTransform.getTranslation().getY()/cm,2));

  if (isEggShapedHyperK){
    geoFile << setw(8)<< 0;
    geoFile << setw(8)<< 0;
  }else{
    geoFile << "Detector radius & height "
	    << setw(8) << innerradius
	    << setw(8)<<WCCylInfo[2];
  }
  geoFile << G4endl;
  
  geoFile << "Type 1 ID PMT number & size "
	  << setw(10)<<totalNumPMTs
	  << setw(8)<<WCPMTSize << setw(4)  <<G4endl;

  geoFile << "Type 2 ID PMT number & size "
	  << setw(10)<<totalNumPMTs2
	  << setw(8)<<WCPMTSize2 << setw(4)  <<G4endl;
  
  geoFile << "OD PMT number & size        "
	  << setw(10)<<totalNumODPMTs
	  << setw(8)<<WCPMTODRadius << setw(4)  <<G4endl;

  geoFile << "Centre offset "
	  << setw(8) << WCOffset(0)
	  << setw(8) << WCOffset(1)
	  << setw(8) << WCOffset(2)
	  << G4endl;
  
  //G4double maxZ=0.0;// used to tell if pmt is on the top/bottom cap
  //G4double minZ=0.0;// or the barrel
  G4int cylLocation = -1;    // Initialize to -1 for error-detection.


  // clear before add new stuff in
  for (unsigned int i=0;i<fpmts.size();i++){
    delete fpmts.at(i);
  }
  fpmts.clear();

  // Grab the tube information from the tubeID Map and dump to file.
  // Reverse loop here, so that cycloc for mPMTs can be made correct
  //  (It is calculated to first order  based on orientation normal to the caps,
  //   which is only correct for the centre 19th 3" PMT in an mPMT module)
  for ( int tubeID = totalNumPMTs; tubeID >= 1; tubeID--){
    G4Transform3D newTransform = tubeIDMap[tubeID];

    // Get tube orientation vector
    G4Vector3D nullOrient = G4Vector3D(0,0,1);
    G4Vector3D pmtOrientation = newTransform * nullOrient;
    // for WCTE type mPMT, center PMT has ID=1
    G4Vector3D pmtOrientation2 = mPMTIDMap[tubeID].second == nID_PMTs ? tubeIDMap[tubeID-nID_PMTs+1] * nullOrient : pmtOrientation ;
    //cyl_location cylLocation = tubeCylLocation[tubeID];

    // Figure out if pmt is on top/bottom or barrel
    // print key: 0-top, 1-barrel, 2-bottom
    if(mPMTIDMap[tubeID].second == 0 || mPMTIDMap[tubeID].second == nID_PMTs) {
      //this only works for non-mPMT (0)
      // or if the highest numbered tube of mPMT is at centre of mPMT module (e.g. on cap points in ±z direction)
      if (!isNuPrism)
      {
        if (pmtOrientation.z()==1.0 || pmtOrientation2.z()==1.0)//bottom
          {cylLocation=2;}
        else if (pmtOrientation.z()==-1.0 || pmtOrientation2.z()==-1.0)//top
          {cylLocation=0;}
        else // barrel
          {cylLocation=1;}
      }
      else // y-axis=vertical axis for nuPrism type detector
      {
        if (pmtOrientation.y()==1.0 || pmtOrientation2.y()==1.0)//bottom
          {cylLocation=2;}
        else if (pmtOrientation.y()==-1.0 || pmtOrientation2.y()==-1.0)//top
          {cylLocation=0;}
        else // barrel
          {cylLocation=1;}
      }
    }

    geoFile.precision(9);
    geoFile << setw(7) << tubeID
	    << " " << setw(5) << mPMTIDMap[tubeID].first
	    << " " << setw(3) << mPMTIDMap[tubeID].second 
 	    << " " << setw(15) << newTransform.getTranslation().getX()/cm
 	    << " " << setw(15) << newTransform.getTranslation().getY()/cm
 	    << " " << setw(15) << newTransform.getTranslation().getZ()/cm
	    << " " << setw(12) << pmtOrientation.x()
	    << " " << setw(12) << pmtOrientation.y()
	    << " " << setw(12) << pmtOrientation.z()
 	    << " " << setw(3) << cylLocation
 	    << G4endl;

     WCSimPmtInfo *new_pmt = new WCSimPmtInfo(cylLocation,
					      newTransform.getTranslation().getX()/cm,
					      newTransform.getTranslation().getY()/cm,
					      newTransform.getTranslation().getZ()/cm,
					      pmtOrientation.x(),
					      pmtOrientation.y(),
					      pmtOrientation.z(),
					      tubeID,
					      mPMTIDMap[tubeID].first,
					      mPMTIDMap[tubeID].second);
     
     fpmts.push_back(new_pmt);

  }
  //reverse the vector of PMTs so that entry 0 = PMT ID 1
  std::reverse(fpmts.begin(), fpmts.end());

  //Record location of the second PMT type for the hybrid configuration
  for (unsigned int i=0;i<fpmts2.size();i++){
    delete fpmts2.at(i);
  }
  fpmts2.clear();

  // Grab the tube information from the ID PMT 2 tubeID Map and dump to file.
  cylLocation = -1;  // Reinitialize to -1 for error-detection.
  for ( int tubeID = totalNumPMTs2; tubeID >= 1; tubeID--){
    G4Transform3D newTransform = tubeIDMap2[tubeID];

    // Get tube orientation vector
    G4Vector3D nullOrient = G4Vector3D(0,0,1);
    G4Vector3D pmtOrientation = newTransform * nullOrient;
    // for WCTE type mPMT, center PMT has ID=1
    G4Vector3D pmtOrientation2 = mPMTIDMap2[tubeID].second == nID_PMTs2 ? tubeIDMap2[tubeID-nID_PMTs2+1] * nullOrient : pmtOrientation ;
    //cyl_location cylLocation = tubeCylLocation[tubeID];

    // Figure out if pmt is on top/bottom or barrel
    // print key: 0-top, 1-barrel, 2-bottom
    if(mPMTIDMap2[tubeID].second == 0 || mPMTIDMap2[tubeID].second == nID_PMTs2) {
      //this only works for non-mPMT (0)
      // or if the highest numbered tube of mPMT is at centre of mPMT module (e.g. on cap points in ±z direction)
      if (pmtOrientation*newTransform.getTranslation() > 0)//veto pmt
	      {cylLocation=10;}
      else if (!isNuPrism)
      {
        if (pmtOrientation.z()==1.0|| pmtOrientation2.z()==1.0)//bottom
          {cylLocation=8;}
        else if (pmtOrientation.z()==-1.0 || pmtOrientation2.z()==-1.0)//top
          {cylLocation=6;}
        else // barrel
          {cylLocation=7;}
      }
      else // y-axis=vertical axis for nuPrism type detector
      {
        if (pmtOrientation.y()==1.0|| pmtOrientation2.y()==1.0)//bottom
          {cylLocation=8;}
        else if (pmtOrientation.y()==-1.0 || pmtOrientation2.y()==-1.0)//top
          {cylLocation=6;}
        else // barrel
          {cylLocation=7;}
      }
    }
    

    geoFile.precision(9);
     geoFile << setw(7) << tubeID
	     << " " << setw(5) << mPMTIDMap2[tubeID].first
	     << " " << setw(3) << mPMTIDMap2[tubeID].second 
 	    << " " << setw(15) << newTransform.getTranslation().getX()/cm
 	    << " " << setw(15) << newTransform.getTranslation().getY()/cm
 	    << " " << setw(15) << newTransform.getTranslation().getZ()/cm
	    << " " << setw(12) << pmtOrientation.x()
	    << " " << setw(12) << pmtOrientation.y()
	    << " " << setw(12) << pmtOrientation.z()
 	    << " " << setw(3) << cylLocation
 	    << G4endl;
     
     WCSimPmtInfo *new_pmt = new WCSimPmtInfo(cylLocation,
					      newTransform.getTranslation().getX()/cm,
					      newTransform.getTranslation().getY()/cm,
					      newTransform.getTranslation().getZ()/cm,
					      pmtOrientation.x(),
					      pmtOrientation.y(),
					      pmtOrientation.z(),
					      tubeID,
					      mPMTIDMap2[tubeID].first,
					      mPMTIDMap2[tubeID].second);
     
     fpmts2.push_back(new_pmt);
  }//loop over PMT type 2
  //reverse the vector of PMTs so that entry 0 = PMT ID 1
  std::reverse(fpmts2.begin(), fpmts2.end());

  // Record locations of OD PMTs to file ffODpmts variables
  for (unsigned int i=0;i<fODpmts.size();i++){
    delete fODpmts.at(i);
  }
  fODpmts.clear();

  // Grab the tube information from the OD tubeID Map and dump to file.
  cylLocation = -1;  // Reinitialize to -1 for error-detection.
  for ( int tubeID = totalNumODPMTs; tubeID >= 1; tubeID--){
    G4Transform3D newTransform = ODtubeIDMap[tubeID];

    // Get tube orientation vector
    G4Vector3D nullOrient = G4Vector3D(0,0,1);
    G4Vector3D pmtOrientation = newTransform * nullOrient;
    //cyl_location cylLocation = tubeCylLocation[tubeID];

    // TODO: make these record something sensible for the OD
    // 3-topOD, 4-barrelOD, 5-bottomOD
    if(mPMTODMap[tubeID].second == 0 || mPMTODMap[tubeID].second == 19) {
      //this only works for non-mPMT (0) or centre tube of mPMT (19)
      if (!isNuPrism)
      {
        if (pmtOrientation.z()==1.0) //top OD
          {cylLocation=5;}
        else if (pmtOrientation.z()==-1.0) //bottom OD
          {cylLocation=3;}
        else // barrel OD
          {cylLocation=4;}
      }
      else // y-axis=vertical axis for nuPrism type detector
      {
        if (pmtOrientation.y()==1.0) //top OD
          {cylLocation=5;}
        else if (pmtOrientation.y()==-1.0) //bottom OD
          {cylLocation=3;}
        else // barrel OD
          {cylLocation=4;}
      }
    }

    geoFile.precision(9);
    geoFile << setw(7) << tubeID
	    << " " << setw(5) << mPMTODMap[tubeID].first
	    << " " << setw(3) << mPMTODMap[tubeID].second
            << " " << setw(15) << newTransform.getTranslation().getX()/CLHEP::cm
            << " " << setw(15) << newTransform.getTranslation().getY()/CLHEP::cm
            << " " << setw(15) << newTransform.getTranslation().getZ()/CLHEP::cm
            << " " << setw(12) << pmtOrientation.x()
            << " " << setw(12) << pmtOrientation.y()
            << " " << setw(12) << pmtOrientation.z()
            << " " << setw(3) << cylLocation
            << G4endl;

    WCSimPmtInfo *new_pmt = new WCSimPmtInfo(cylLocation,
                                             newTransform.getTranslation().getX()/CLHEP::cm,
                                             newTransform.getTranslation().getY()/CLHEP::cm,
                                             newTransform.getTranslation().getZ()/CLHEP::cm,
                                             pmtOrientation.x(),
                                             pmtOrientation.y(),
                                             pmtOrientation.z(),
                                             tubeID,
					     mPMTODMap[tubeID].first,
					     mPMTODMap[tubeID].second);

    fODpmts.push_back(new_pmt);

  }//loop over OD PMTs
  //reverse the vector of PMTs so that entry 0 = PMT ID 1
  std::reverse(fODpmts.begin(), fODpmts.end());

  geoFile.close();

  G4cout << "Geofile written" << G4endl;
} 

// Code for traversing the geometry tree.  This code is very general you pass
// it a function and it will call the function with the information on each
// object it finds.
//
// The traversal code comes from a combination of me/G4Lab project &
// from source/visualization/modeling/src/G4PhysicalVolumeModel.cc
//
// If you are trying to understand how passing the function works you need
// to understand pointers to member functions...
//
// Also notice that DescriptionFcnPtr is a (complicated) typedef.
//

void WCSimDetectorConstruction::TraverseReplicas
(G4VPhysicalVolume* aPV, int aDepth, const G4Transform3D& aTransform,
 DescriptionFcnPtr registrationRoutine)
{
  // Recursively visit all of the geometry below the physical volume
  // pointed to by aPV including replicas.
  
  G4ThreeVector     originalTranslation = aPV->GetTranslation();
  G4RotationMatrix* pOriginalRotation   = aPV->GetRotation();
  
  if (aPV->IsReplicated() ) 
  {
    EAxis    axis;
    G4int    nReplicas;
    G4double width, offset;
    G4bool   consuming;

    aPV->GetReplicationData(axis,nReplicas,width,offset,consuming);

    for (int n = 0; n < nReplicas; n++) 
      {
	//TF: for mPMT parameterization
	// because kZAxis is a dummy value altered by Parameterised volume in MultiPMTParam
	if(aPV->IsParameterised() ){
	  //G4cout << "Replica No. " << n << G4endl;
	  //use the actual translation and rotation of the replica.
	  (aPV->GetParameterisation())->ComputeTransformation(n,aPV); 
	}
	
	else{

	  switch(axis) {
	  default:
	  case kXAxis:
	    aPV->SetTranslation(G4ThreeVector
				(-width*(nReplicas-1)*0.5+n*width,0,0));
	    aPV->SetRotation(0);
	    break;
	  case kYAxis:
	    aPV->SetTranslation(G4ThreeVector
				(0,-width*(nReplicas-1)*0.5+n*width,0));
	    aPV->SetRotation(0);
	    break;
	  case kZAxis:
	    aPV->SetTranslation(G4ThreeVector
				(0,0,-width*(nReplicas-1)*0.5+n*width));
	    aPV->SetRotation(0);
	    break;
	  case kRho:
	    //Lib::Out::putL("GeometryVisitor::visit: WARNING:");
	    //Lib::Out::putL(" built-in replicated volumes replicated");
	    //Lib::Out::putL(" in radius are not yet properly visualizable.");
	    aPV->SetTranslation(G4ThreeVector(0,0,0));
	    aPV->SetRotation(0);
	    break;
	  case kPhi:
	    {
	      G4RotationMatrix rotation;
	      rotation.rotateZ(-(offset+(n+0.5)*width));
	      // Minus Sign because for the physical volume we need the
	      // coordinate system rotation.
	      aPV->SetTranslation(G4ThreeVector(0,0,0));
	      aPV->SetRotation(&rotation);
	    }
	    break;
	    
	  } // axis switch
	}
      DescribeAndDescendGeometry(aPV, aDepth, n, aTransform, 
				 registrationRoutine);

    }   // num replicas for loop
  }     // if replicated
  else
    DescribeAndDescendGeometry(aPV, aDepth, aPV->GetCopyNo(), aTransform, 
			       registrationRoutine);
  
  // Restore original transformation...
  aPV->SetTranslation(originalTranslation);
  aPV->SetRotation(pOriginalRotation);
}

void WCSimDetectorConstruction::DescribeAndDescendGeometry
(G4VPhysicalVolume* aPV ,int aDepth, int replicaNo, 
 const G4Transform3D& aTransform,  DescriptionFcnPtr registrationRoutine)
{

  // Calculate the new transform relative to the old transform
  G4Transform3D* transform = 
    new G4Transform3D(*(aPV->GetObjectRotation()), aPV->GetTranslation());

  G4Transform3D newTransform = aTransform * (*transform);
  delete transform; 

  /*
  G4cout << aPV->GetObjectRotation()->getPhi() << " " << aPV->GetObjectRotation()->getTheta() << " " 
	    << aPV->GetObjectRotation()->getPsi() << G4endl;
  G4cout << aPV->GetTranslation().x() << " " << aPV->GetTranslation().y() << " " 
	    << aPV->GetTranslation().z() << G4endl;
  */

  // Call the routine we use to print out geometry descriptions, make
  // tables, etc.  The routine was passed here as a parameter.  It needs to
  // be a member function of the class

  (this->*registrationRoutine)(aPV, aDepth, replicaNo, newTransform);

  int nDaughters = aPV->GetLogicalVolume()->GetNoDaughters();
    
  for (int iDaughter = 0; iDaughter < nDaughters; iDaughter++) 
    TraverseReplicas(aPV->GetLogicalVolume()->GetDaughter(iDaughter),
		     aDepth+1, newTransform, registrationRoutine);
}



G4double WCSimDetectorConstruction::GetGeo_Dm(G4int i){
  if (i>=0&&i<=2){
    return WCCylInfo[i];
  }else if(i==3){
    return innerradius;
  }else{
    return 0;
  }
}

// Read PMT positions from file
void WCSimDetectorConstruction::ReadGeometryTableFromFile(){
  pmtPos.clear();
  pmtDir.clear();
  pmtUse.clear();
  pmtType.clear();
  pmtRotaton.clear();
  pmtSection.clear();
  pmtmPMTId.clear();
  nPMTsRead = 0;
  if (readFromTable) ReadGeometryTableFromFile(pmtPositionFile);
  if (readODFromTable) ReadGeometryTableFromFile(odpmtPositionFile);
}

void WCSimDetectorConstruction::ReadGeometryTableFromFile(std::string fname)
{
  std::ifstream Data(fname.c_str(),std::ios_base::in);
  if (!Data)
  {
    G4cout<<"PMT data file "<<fname<<" could not be opened --> Exiting..."<<G4endl;
    exit(-1);
  }
  else
    G4cout<<"PMT data file "<<fname<<" is opened to read positions"<<G4endl;

  std::string str, tmp;
	G4int Column=0;
	while (std::getline(Data, str)) {
		if (str=="#DATASTART") break;
	}
	std::ifstream::pos_type SavePoint = Data.tellg();
	std::getline(Data, str);
	std::istringstream stream(str);
	while (std::getline(stream,tmp,' ')) Column++;
	if (Column!=11)
  {
    G4cerr<<"Number of column = "<<Column<<" which is not equal to 11. "<<G4endl;
    G4cerr<<"Inappropriate input --> Exiting..."<<G4endl;
    exit(-1);
  }
  Data.seekg(SavePoint);

  G4ThreeVector pos, dir;
  double angle;
  int ptype, UsePMT, sect, mpmtid;
  while (!Data.eof())
  {
    Data>>ptype>>mpmtid>>sect>>pos[0]>>pos[1]>>pos[2]>>dir[0]>>dir[1]>>dir[2]>>angle>>UsePMT;
    pmtPos.push_back(pos);
    pmtDir.push_back(dir);
    pmtUse.push_back(UsePMT);
    pmtType.push_back(ptype);
    pmtSection.push_back(sect);
    pmtmPMTId.push_back(mpmtid);
    pmtRotaton.push_back(angle*deg);
    nPMTsRead++;
  }
  Data.close();
}