
#include "vtkObjectFactory.h"
#include "vtkVersion.h"

        

class vtkCustomObjectFactory: public vtkObjectFactory
{
public:  
    
    vtkTypeMacro( vtkCustomObjectFactory, vtkObjectFactory );
    static vtkCustomObjectFactory* New();
    
    virtual vtkObject* CreateObject( const char* classname );
    
    char* GetVTKSourceVersion() { return VTK_VERSION; }
    
    char* GetDescription() { return "Overrides for META project"; }
};


