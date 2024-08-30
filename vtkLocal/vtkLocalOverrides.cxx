#include "vtkJCNArrayCalculator.h"
#include "vtkJCNMergeFields.h"
#include "vtkJCNMergePoints.h"
#include "vtkJCNSplitField.h"

#include "vtkObjectFactory.h"
#include "vtkVersion.h"

class vtkCustomObjectFactory: public vtkObjectFactory
{
public:  
    
    vtkTypeMacro( vtkCustomObjectFactory, vtkObjectFactory );
    static vtkCustomObjectFactory* New();
    
    virtual vtkObject* CreateObject( const char* classname )
    {

        #define OVERRIDE( oldcn, newcn ) \
        if( !strcmp( classname, #oldcn ) ) { \
            return newcn::New(); \
        }

        OVERRIDE( vtkMergeFields, vtkJCNMergeFields );
        OVERRIDE( vtkArrayCalculator, vtkJCNArrayCalculator );
        OVERRIDE( vtkMergePoints, vtkJCNMergePoints );
        OVERRIDE( vtkSplitField, vtkJCNSplitField );
        
        return NULL;
    }
    
    char* GetVTKSourceVersion() {
        return VTK_VERSION;
    }
    
    char* GetDescription() {
        return "foo";
    }
};

vtkStandardNewMacro(vtkCustomObjectFactory);

// --------------------------------

__attribute__((constructor)) static void setup() 
{
    std::cerr << "vtkLocal: class overrides enabled\n";

    vtkCustomObjectFactory* cof = vtkCustomObjectFactory::New();
    vtkObjectFactory::RegisterFactory( cof );
}
