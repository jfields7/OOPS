#ifndef FIELD_MAP_H
#define FIELD_MAP_H
#include <map>
#include <string>
#include <memory>

class ODEData;
class SolverData;
class Grid;
class FieldInfo;

class FieldMap{
  private:
    /**
     * A constant reference to the Grid object this should be sorted by.
     */
    const Grid& mGrid;

    /**
     * A string-sorted map of ODEData objects.
     */
    std::map<std::string, std::shared_ptr<ODEData>> fields;
    /**
     * A string-sorted map of SolverData objects.
     */
    std::map<std::string, std::shared_ptr<SolverData>> solverFields;
  public:
    /**
     * The FieldMap constructor. The Grid object must be specified at
     * construction and is immutable.
     */
    FieldMap(const Grid& grid);
    /**
     * An alternate FieldMap constructor. It takes a map of FieldInfo
     * objects and automatically adds them to the FieldMap.
     */
    FieldMap(const Grid& grid, std::map<std::string, FieldInfo>& fields);

    /**
     * The FieldMap destructor. It clears all the fields.
     */
     ~FieldMap();

    /**
     * Add a new field to the map. Fields that will be evolved inside a
     * Solver should include the number of stages in the solver. Setting
     * stages to a nonzero number indicates that the resultant field is
     * going to be evolved by a Solver object and should be a SolverData
     * object rather than an ODEData object.
     * @param name The unique name for this field.
     * @param nEqs The number of equations in the new field.
     * @param stages An optional parameter (0 by default) describing the
     *        number of stages the SolverData object needs.
     */
    void addField(std::string name, unsigned int nEqs, unsigned int stages);

    /**
     * Remove a field from the map.
     * @param The name of the field.
     * @exception Throws std::out_of_range when a field does not
     *            exist.
     */
    void removeField(std::string name);

    /**
     * Check if a particular field exists.
     * @param The name of the field.
     * @return True if the field exists, 
     *         false otherwise.
     */
    bool hasField(std::string name);

    /**
     * Get the number of different fields.
     * @return The total number of fields stored in the map.
     */
    inline unsigned int getNumFields(){
      return fields.size();
    }

    /**
     * Check if a particular field is a SolverData field or just an
     * ODEData field.
     * @param The name of the field.
     * @return True if the field is a SolverData object, false if
     *         otherwise.
     * @exception Throws std::out_of_range when a field does not
     *            exist.
     */
    bool isSolverField(std::string name);

    inline const Grid& getGrid() const{
      return mGrid;
    }

    /**
     * Return a field as a SolverData field.
     * @param The name of the field
     * @return A reference to a shared_ptr cast to a SolverData object.
     * @exception Throws std::out_of_range when a field does not exist.
     */
    inline std::shared_ptr<SolverData>& getSolverField(std::string name){
      return solverFields.at(name);
    }

    /**
     * Overload the [] operator for easy access.
     * @param name The name of the field to access.
     * @return A reference to the ODEData field.
     * @exception Throws std::out_of_range when 
     *            a field does not exist.
     */
    inline std::shared_ptr<ODEData>& operator[](std::string name){
      return fields.at(name);
    }

    /**
     * Overload the < operator so that fieldMap objects
     * can be compared by their Grid objects. This enables
     * them to be spatially sorted.
     */
    bool operator < (const FieldMap& other);
};

#endif
