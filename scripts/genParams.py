# Read in a JSON setup file and generate a parameter file for it.
import json
import sys

def addHeader(f, name):
  f.write("#ifndef %s_PARAMETERS_H\n" % (name.upper()));
  f.write("#define %s_PARAMETERS_H\n\n" % (name.upper()));
  f.write("#include <parameters.h>\n\n");
  f.write("class %sParameters : public Parameters {\n" % name);

def addEnum(f, name, params):
  f.write("    enum %s{\n" % name);
  for s in params:
    f.write("      %s,\n" % s.upper());
  f.write("    };\n\n");

def addFooter(f):
  f.write("};\n\n");
  f.write("#endif\n");

def addGetter(f, name, vtype):
  f.write("    inline %s get%s(){\n" % (vtype,name));
  f.write("      return m%s;\n" % name);
  f.write("    }\n\n");

def addEnumSetter(f, name):
  f.write("    inline void set%s(%s val){\n" % (name, name));
  f.write("      m%s = val;\n" % (name));
  f.write("    }\n\n");

def addSetter(f, name, vtype):
  f.write("    inline void set%s(%s %s){\n" % (name, vtype, name));
  f.write("      m%s = %s;\n" % (name, name));
  f.write("    }\n\n");

def addVariable(f, name, vtype):
  f.write("    %s m%s;\n" % (vtype, name));

def addPublic(f):
  f.write("  public:\n");

def addPrivate(f):
  f.write("  private:\n");

def addConstructor(f, obj):
  f.write("    %sParameters() : Parameters(%i){\n" % (obj["name"], obj["id"]));

  members = obj["members"];
  for var in members:
    f.write("      m" + var["name"] + " = " + str(var["default"]) + ";\n");

  f.write("    }\n\n");

def addEnums(f, obj):
  members = obj["members"];
  for var in members:
    if var["type"] == "enum":
      addEnum(f, var["name"], var["value"]);

def writeClass(obj):
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

def writeParser(obj):
  f = open(obj["name"].lower() + "parser.h","w");

  # Add the header.
  f.write("#ifndef %s_PARSER_H\n" % obj["name"].upper());
  f.write("#define %s_PARSER_H\n\n" % obj["name"].upper());
  f.write("#include <paramparser.h>\n");
  f.write("#include <%sparameters.h>\n\n" % obj["name"].lower());

  # Start writing the class.
  f.write("class %sParser : ParamParser{\n" % obj["name"]);
  addPublic(f);
  f.write("    %sParser() : ParamParser(%i){}\n\n" % (obj["name"], obj["id"]));
  f.write("    virtual ~%sParser(){}\n\n" % obj["name"]);

  # Write the getParameters() function.
  f.write("    virtual Parameters& getParameters(std::string fname){\n");
  f.write("      %sParameters params = %sParameters();\n" % (obj["name"], obj["name"]));
  f.write("      updateParameters(fname, &params);\n");
  f.write("      return params;\n");
  f.write("    }\n\n");

  f.write("    virtual void updateParameters(std::string fname, Parameters* params);\n");

  addFooter(f);
  f.close();

  # Write the .cpp file.
  f = open(obj["name"].lower() + "parser.cpp","w");

  f.write("#include <%sparser.h>\n" % obj["name"].lower());
  f.write("#include <iostream>\n\n");
  f.write("using namespace std;\n\n");

  # Write the updateParameters() function.
  f.write("void %sParser::updateParameters(string fname, Parameters* params){\n" % obj["name"]);
  f.write("  if(!checkId(params)){\n");
  f.write("    return;\n");
  f.write("  }\n");
  f.write("  reader.clearData();\n");
  f.write("  ParamReader::ParamResult result = reader.readFile(fname);\n");
  f.write("  if(result != ParamReader::SUCCESS){\n");
  f.write("    cout << \"An error occurred while trying to read \" << fname << \".\\n\";\n");
  f.write("    return;\n");
  f.write("  }\n");
  f.write("  %sParameters *pars = (%sParameters*) params;\n\n" % (obj["name"], obj["name"]));
  f.write("  if(!reader.hasSection(string(\"%s\"))){\n" % obj["name"]);
  f.write("    return;\n");
  f.write("  }\n\n");

  # Loop through all the different parameters and check their types.
  for var in obj["members"]:
    f.write("  if(reader.hasParameter(string(\"%s\"),string(\"%s\"))){\n" % (obj["name"], var["name"]));
    # Parse a double.
    if var["type"] == "double":
      f.write("    double result = reader.readAsDouble(string(\"%s\"),string(\"%s\"));\n" % (obj["name"], var["name"]));
      f.write("    if(result <= %e && result >= %e){\n" % (var["max"], var["min"]));
      f.write("      pars->set%s(result);\n" % var["name"]);
      f.write("    }\n");
      f.write("    else{\n");
      f.write("      cout << \"Parameter %s out of range.\\n\";\n");
      f.write("    }\n");
      f.write("  }\n\n");
    # Parse an int.
    elif var["type"] == "int":
      f.write("    int result = reader.readAsInt(string(\"%s\"),string(\"%s\"));\n" % (obj["name"], var["name"]));
      f.write("    if(result <= %d && result >= %d){\n" % (var["max"], var["min"]));
      f.write("      pars->set%s(result);\n" % var["name"]);
      f.write("    }\n");
      f.write("    else{\n");
      f.write("      cout << \"Parameter %s out of range.\\n\";\n");
      f.write("    }\n");
      f.write("  }\n\n");
    # Parse a string.
    elif var["type"] == "string":
      f.write("    string result = reader.readAsString(string(\"%s\"),string(\"%s\"));\n" % (obj["name"], var["name"]));
      f.write("    if(result.compare(\"NULL\") != 0){\n");
      f.write("      pars->set%s(result);\n" % var["name"]);
      f.write("    }\n");
      f.write("    else{\n");
      f.write("      cout << \"Parameter %s is not a string.\\n\";\n");
      f.write("    }\n");
      f.write("  }\n\n");
    # Parse an enumerator.
    elif var["type"] == "enum":
      f.write("    string result = reader.readAsString(string(\"%s\"),string(\"%s\"));\n" % (obj["name"], var["name"]));
      f.write("    "); # The else if construction means we need to add the space separately.
      for s in var["value"]:
        f.write("if(result.compare(\"%s\") == 0){\n" % s);
        f.write("      pars->set%s(%sParameters::%s);\n" % (var["name"], obj["name"], s));
        f.write("    }\n");
        f.write("    else ");
      f.write("{\n");
      f.write("      cout << \"The value for parameter %s is out of range.\\n\";\n" % var["name"]);
      f.write("    }\n");
      f.write("  }\n\n");

  f.write("}\n\n");
  f.close();


# Main section of the script.
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
      writeClass(data[key]);
      writeParser(data[key]);

    
  except json.JSONDecodeError as err:
    print("There was an error decoding " + sys.argv[i] + ":");
    print(" " + err.msg + " at line " + str(err.lineno))
  except:
    print("Could not read file " + sys.argv[i] + ".");
