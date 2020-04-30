#include <paramreader.h>
#include <cstdio>
#include <string>

using namespace std;

void printTest(char *name, bool success){
  if(success){
    printf("\033[1;32m%s test passed.\033[0m\n\n",name);
  }
  else{
    printf("\033[1;31m%s test failed.\033[0m\n\n",name);
  }
}

int main(int argc, char *argv[]){
  ParamReader reader;

  ParamReader::ParamResult result = reader.readFile("test.ini");

  // Now we validate the results.
  printTest("Load File", result == ParamReader::SUCCESS);

  // Make sure the correct sections exist.
  bool r = reader.hasSection("GLOBAL");
  r = reader.hasSection("NEXT SECTION") && r;
  r = !reader.hasSection("BOB") && r;
  printTest("Check Sections", r);

  // Make sure we have the right parameters.
  r = reader.hasParameter("GLOBAL","bob");
  r = reader.hasParameter("GLOBAL","steve") && r;
  r = reader.hasParameter("GLOBAL","old") && r;
  r = reader.hasParameter("NEXT SECTION","foo") && r;
  printTest("Loaded Parameters 1", r);
  // Make sure the parameters aren't in the wrong place.
  r = !reader.hasParameter("NEXT SECTION","bob");
  r = !reader.hasParameter("NEXT SECTION","steve") && r;
  r = !reader.hasParameter("NEXT SECTION","old") && r;
  r = !reader.hasParameter("GLOBAL","foo") && r;
  printTest("Loaded Parameters 2", r);

  // Now check the parameter values.
  int i = reader.readAsInt("GLOBAL","bob");
  printTest("Read Integer", i == 10);
  i = reader.readAsInt("GLOBAL","old");
  printTest("Read Integer 2", i == 1);
  string str = reader.readAsString("GLOBAL","steve");
  printTest("Read String",str.compare("derp") == 0);
  double d = reader.readAsDouble("NEXT SECTION","foo");
  printTest("Read Double", d == 0.1);
  // Make sure that the parameters can't be read as something else.
  i = reader.readAsInt("NEXT SECTION","foo");
  i += reader.readAsInt("GLOBAL","steve");
  printTest("Read Bad Integer",i == 0);

  // Clear everything.
  reader.clearData();

  // Unsectioned parameter test
  result = reader.readFile("unsectioned.ini");
  printTest("Unsectioned Parameter",result == ParamReader::UNSECTIONED_PARAMETER);

  reader.clearData();

  // Multiple definition test
  result = reader.readFile("multiple.ini");
  printTest("Multiple Definition", result == ParamReader::MULTIPLE_DEFINITIONS);


  

  return 0;
}
