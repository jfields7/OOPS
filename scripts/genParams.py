# Read in a JSON setup file and generate a parameter file for it.
import json
import sys

def addHeader(f, name):
  f.write("#ifndef %s_PARAMETERS_H\n" % (name.upper()));
  f.write("#define %s_PARAMETERS_H\n\n" % (name.upper()));
  f.write("#include <parameters.h>\n\n");
  f.write("class %sParameters : public Parameters {\n" % name);

def addEnum(f, name, params):
  f.write("   enum %s{\n" % name);
  for s in params:
    f.write("     %s,\n" % s.upper());
  f.write("   };\n\n");

def addFooter(f):
  f.write("};\n\n");
  f.write("#endif\n");

def addGetter(f, name, vtype):
  f.write("   inline %s get%s(){\n" % (vtype,name));
  f.write("     return m%s;\n" % name);
  f.write("   }\n\n");

def addEnumSetter(f, name):
  f.write("   inline void set%s(%s val){\n" % (name, name));
  f.write("     m%s = val;\n" % (name));
  f.write("   }\n\n");

def addSetter(f, name, vtype):
  f.write("   inline void set%s(%s %s){\n" % (name, vtype, name));
  f.write("     m%s = %s;\n" % (name, name));
  f.write("   }\n\n");

def addVariable(f, name, vtype):
  f.write("   %s m%s;\n" % (vtype, name));

def addPublic(f):
  f.write(" public:\n");

def addPrivate(f):
  f.write(" private:\n");

def addConstructor(f, obj):
  f.write("   %sParameters() : Parameters(%i){\n" % (obj["name"], obj["id"]));

  members = obj["members"];
  for var in members:
    f.write("     m" + var["name"] + " = " + str(var["default"]) + ";\n");

  f.write("   }\n\n");

def addEnums(f, obj):
  members = obj["members"];
  for var in members:
    if var["type"] == "enum":
      addEnum(f, var["name"], var["value"]);

def parseClass(obj):
  f = open(obj["name"].lower() + "parameters.h","w");
  addHeader(f,obj["name"]);

  addPublic(f);

  addEnums(f,obj);
  
  addConstructor(f, obj);

  members = obj["members"];
  for var in members:
    if var["type"] == "enum":
      addEnumSetter(f, var["name"]);
      addGetter(f, var["name"], var["name"]);
    else:
      addSetter(f, var["name"], var["type"]);
      addGetter(f, var["name"], var["type"]);

  addPrivate(f);
  for var in members:
    if var["type"] == "enum":
      addVariable(f, var["name"], var["name"]);
    else:
      addVariable(f, var["name"], var["type"]);

  addFooter(f);
  f.close();

if len(sys.argv) == 1:
  print("Please call genParams.py with a filename or list of filenames.");
  quit();

for i in range(1,len(sys.argv)):
  try:
    f = open(sys.argv[i],"r");
    data = json.load(f);
    f.close();

    # We assume that the first level of objects is always a class name.
    for key in data:
      parseClass(data[key]);

    
  except json.JSONDecodeError as err:
    print("There was an error decoding " + sys.argv[i] + ":");
    print(" " + err.msg + " at line " + str(err.lineno))
  except:
    print("Could not read file " + sys.argv[i] + ".");
