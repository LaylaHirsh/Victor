#include <Group.h>
#include <vector3.h>
#include <XyzLoader.h>
#include <XyzSaver.h>
#include <vglMath.h>

using namespace Biopool;

vg_ieee64 degreesToRadian(const vg_ieee64& deg)
{
return (deg / 180.0) * M_PI;
}

int main()
{ 
  cout << "Start" << endl;

  Group g, g2;
  ifstream inFile("test-grp.xyz");
  if (!inFile)
    ERROR("File not found.", exception);

  XyzLoader il(inFile);
  g.load(il);

  ifstream inFile2("test-grp.xyz");
  if (!inFile2)
    ERROR("File not found.", exception);
  XyzLoader il2(inFile2);
  g2.load(il2);

  g2.setType("2-C3H8");

  cout << "-------------------------------------------------------------\n";

  XyzSaver is(cout, 0);
  g.save(is);

  cout << "-------------------------------------------------------------\n";
  
  g2.save(is);
  g2.bindIn(g2[0], g, g[6]);

  cout << "-------------------------------------------------------------\n";
  g2.save(is);
  
  vgVector3<double> tmp(0,0,1);
  vgMatrix3<double> rotationMatrix = 
          vgMatrix3<double>::createRotationMatrix(tmp, degreesToRadian(180));
  g2[0].addRot(rotationMatrix);

  cout << "-------------------------------------------------------------\n";
  g2.save(is);

  g2.addRot(rotationMatrix);

  cout << "-------------------------------------------------------------\n";
  g2.save(is);

  cout << "-------------------------------------------------------------\n";

  cout << "superior: \t g.a1= " << g[1].getSuperior().getType() 
       << "\t g2.a3= " << g2[3].getSuperior().getType() << endl;

  return 0;
}
