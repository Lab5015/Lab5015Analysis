// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME dOdIlibdImydict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "./interface/SetTDRStyle.h"
#include "./interface/AnalysisUtils.h"
#include "./interface/FitUtils.h"
#include "./interface/Na22SpectrumAnalyzer.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_EventClass(void *p = 0);
   static void *newArray_EventClass(Long_t size, void *p);
   static void delete_EventClass(void *p);
   static void deleteArray_EventClass(void *p);
   static void destruct_EventClass(void *p);
   static void streamer_EventClass(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::EventClass*)
   {
      ::EventClass *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::EventClass >(0);
      static ::ROOT::TGenericClassInfo 
         instance("EventClass", ::EventClass::Class_Version(), "AnalysisUtils.h", 19,
                  typeid(::EventClass), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::EventClass::Dictionary, isa_proxy, 16,
                  sizeof(::EventClass) );
      instance.SetNew(&new_EventClass);
      instance.SetNewArray(&newArray_EventClass);
      instance.SetDelete(&delete_EventClass);
      instance.SetDeleteArray(&deleteArray_EventClass);
      instance.SetDestructor(&destruct_EventClass);
      instance.SetStreamerFunc(&streamer_EventClass);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::EventClass*)
   {
      return GenerateInitInstanceLocal((::EventClass*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::EventClass*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_ModuleEventClass(void *p = 0);
   static void *newArray_ModuleEventClass(Long_t size, void *p);
   static void delete_ModuleEventClass(void *p);
   static void deleteArray_ModuleEventClass(void *p);
   static void destruct_ModuleEventClass(void *p);
   static void streamer_ModuleEventClass(TBuffer &buf, void *obj);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::ModuleEventClass*)
   {
      ::ModuleEventClass *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::ModuleEventClass >(0);
      static ::ROOT::TGenericClassInfo 
         instance("ModuleEventClass", ::ModuleEventClass::Class_Version(), "AnalysisUtils.h", 62,
                  typeid(::ModuleEventClass), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::ModuleEventClass::Dictionary, isa_proxy, 16,
                  sizeof(::ModuleEventClass) );
      instance.SetNew(&new_ModuleEventClass);
      instance.SetNewArray(&newArray_ModuleEventClass);
      instance.SetDelete(&delete_ModuleEventClass);
      instance.SetDeleteArray(&deleteArray_ModuleEventClass);
      instance.SetDestructor(&destruct_ModuleEventClass);
      instance.SetStreamerFunc(&streamer_ModuleEventClass);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::ModuleEventClass*)
   {
      return GenerateInitInstanceLocal((::ModuleEventClass*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::ModuleEventClass*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr EventClass::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *EventClass::Class_Name()
{
   return "EventClass";
}

//______________________________________________________________________________
const char *EventClass::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EventClass*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int EventClass::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::EventClass*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *EventClass::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EventClass*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *EventClass::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::EventClass*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr ModuleEventClass::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *ModuleEventClass::Class_Name()
{
   return "ModuleEventClass";
}

//______________________________________________________________________________
const char *ModuleEventClass::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ModuleEventClass*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int ModuleEventClass::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::ModuleEventClass*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *ModuleEventClass::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ModuleEventClass*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *ModuleEventClass::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::ModuleEventClass*)0x0)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
void EventClass::Streamer(TBuffer &R__b)
{
   // Stream an object of class EventClass.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      { TString R__str; R__str.Streamer(R__b); stepLabel = R__str.Data(); }
      { TString R__str; R__str.Streamer(R__b); ch1 = R__str.Data(); }
      { TString R__str; R__str.Streamer(R__b); ch2 = R__str.Data(); }
      { TString R__str; R__str.Streamer(R__b); label1 = R__str.Data(); }
      { TString R__str; R__str.Streamer(R__b); label2 = R__str.Data(); }
      { TString R__str; R__str.Streamer(R__b); label12 = R__str.Data(); }
      R__b >> x;
      R__b >> y;
      R__b >> isBar1;
      R__b >> isBar2;
      R__b >> isBarSide1;
      R__b >> isBarSide2;
      R__b >> isHorizontal1;
      R__b >> isHorizontal2;
      R__b >> qfine1;
      R__b >> qfine1L;
      R__b >> qfine1R;
      R__b >> qfine2;
      R__b >> qfine2L;
      R__b >> qfine2R;
      R__b >> tot1;
      R__b >> tot1L;
      R__b >> tot1R;
      R__b >> tot2;
      R__b >> tot2L;
      R__b >> tot2R;
      R__b >> energy1;
      R__b >> energy1L;
      R__b >> energy1R;
      R__b >> energy2;
      R__b >> energy2L;
      R__b >> energy2R;
      R__b >> time1;
      R__b >> time2;
      R__b.CheckByteCount(R__s, R__c, EventClass::IsA());
   } else {
      R__c = R__b.WriteVersion(EventClass::IsA(), kTRUE);
      TObject::Streamer(R__b);
      { TString R__str = stepLabel.c_str(); R__str.Streamer(R__b);}
      { TString R__str = ch1.c_str(); R__str.Streamer(R__b);}
      { TString R__str = ch2.c_str(); R__str.Streamer(R__b);}
      { TString R__str = label1.c_str(); R__str.Streamer(R__b);}
      { TString R__str = label2.c_str(); R__str.Streamer(R__b);}
      { TString R__str = label12.c_str(); R__str.Streamer(R__b);}
      R__b << x;
      R__b << y;
      R__b << isBar1;
      R__b << isBar2;
      R__b << isBarSide1;
      R__b << isBarSide2;
      R__b << isHorizontal1;
      R__b << isHorizontal2;
      R__b << qfine1;
      R__b << qfine1L;
      R__b << qfine1R;
      R__b << qfine2;
      R__b << qfine2L;
      R__b << qfine2R;
      R__b << tot1;
      R__b << tot1L;
      R__b << tot1R;
      R__b << tot2;
      R__b << tot2L;
      R__b << tot2R;
      R__b << energy1;
      R__b << energy1L;
      R__b << energy1R;
      R__b << energy2;
      R__b << energy2L;
      R__b << energy2R;
      R__b << time1;
      R__b << time2;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_EventClass(void *p) {
      return  p ? new(p) ::EventClass : new ::EventClass;
   }
   static void *newArray_EventClass(Long_t nElements, void *p) {
      return p ? new(p) ::EventClass[nElements] : new ::EventClass[nElements];
   }
   // Wrapper around operator delete
   static void delete_EventClass(void *p) {
      delete ((::EventClass*)p);
   }
   static void deleteArray_EventClass(void *p) {
      delete [] ((::EventClass*)p);
   }
   static void destruct_EventClass(void *p) {
      typedef ::EventClass current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_EventClass(TBuffer &buf, void *obj) {
      ((::EventClass*)obj)->::EventClass::Streamer(buf);
   }
} // end of namespace ROOT for class ::EventClass

//______________________________________________________________________________
void ModuleEventClass::Streamer(TBuffer &R__b)
{
   // Stream an object of class ModuleEventClass.

   UInt_t R__s, R__c;
   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(&R__s, &R__c); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> barID;
      R__b >> Vov;
      R__b >> vth1;
      R__b >> energyL;
      R__b >> energyR;
      R__b >> timeL;
      R__b >> timeR;
      R__b.CheckByteCount(R__s, R__c, ModuleEventClass::IsA());
   } else {
      R__c = R__b.WriteVersion(ModuleEventClass::IsA(), kTRUE);
      TObject::Streamer(R__b);
      R__b << barID;
      R__b << Vov;
      R__b << vth1;
      R__b << energyL;
      R__b << energyR;
      R__b << timeL;
      R__b << timeR;
      R__b.SetByteCount(R__c, kTRUE);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_ModuleEventClass(void *p) {
      return  p ? new(p) ::ModuleEventClass : new ::ModuleEventClass;
   }
   static void *newArray_ModuleEventClass(Long_t nElements, void *p) {
      return p ? new(p) ::ModuleEventClass[nElements] : new ::ModuleEventClass[nElements];
   }
   // Wrapper around operator delete
   static void delete_ModuleEventClass(void *p) {
      delete ((::ModuleEventClass*)p);
   }
   static void deleteArray_ModuleEventClass(void *p) {
      delete [] ((::ModuleEventClass*)p);
   }
   static void destruct_ModuleEventClass(void *p) {
      typedef ::ModuleEventClass current_t;
      ((current_t*)p)->~current_t();
   }
   // Wrapper around a custom streamer member function.
   static void streamer_ModuleEventClass(TBuffer &buf, void *obj) {
      ((::ModuleEventClass*)obj)->::ModuleEventClass::Streamer(buf);
   }
} // end of namespace ROOT for class ::ModuleEventClass

namespace ROOT {
   static TClass *vectorlEvectorlEvectorlEintgRsPgRsPgR_Dictionary();
   static void vectorlEvectorlEvectorlEintgRsPgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEvectorlEintgRsPgRsPgR(void *p = 0);
   static void *newArray_vectorlEvectorlEvectorlEintgRsPgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEvectorlEintgRsPgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEvectorlEintgRsPgRsPgR(void *p);
   static void destruct_vectorlEvectorlEvectorlEintgRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<vector<int> > >*)
   {
      vector<vector<vector<int> > > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<vector<int> > >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<vector<int> > >", -2, "vector", 210,
                  typeid(vector<vector<vector<int> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEvectorlEintgRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<vector<int> > >) );
      instance.SetNew(&new_vectorlEvectorlEvectorlEintgRsPgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEvectorlEintgRsPgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEvectorlEintgRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEvectorlEintgRsPgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEvectorlEintgRsPgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<vector<int> > > >()));

      ::ROOT::AddClassAlternate("vector<vector<vector<int> > >","std::vector<std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >, std::allocator<std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<vector<int> > >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEvectorlEintgRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<vector<int> > >*)0x0)->GetClass();
      vectorlEvectorlEvectorlEintgRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEvectorlEintgRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEvectorlEintgRsPgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<vector<int> > > : new vector<vector<vector<int> > >;
   }
   static void *newArray_vectorlEvectorlEvectorlEintgRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<vector<int> > >[nElements] : new vector<vector<vector<int> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEvectorlEintgRsPgRsPgR(void *p) {
      delete ((vector<vector<vector<int> > >*)p);
   }
   static void deleteArray_vectorlEvectorlEvectorlEintgRsPgRsPgR(void *p) {
      delete [] ((vector<vector<vector<int> > >*)p);
   }
   static void destruct_vectorlEvectorlEvectorlEintgRsPgRsPgR(void *p) {
      typedef vector<vector<vector<int> > > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<vector<int> > >

namespace ROOT {
   static TClass *vectorlEvectorlEvectorlEfloatgRsPgRsPgR_Dictionary();
   static void vectorlEvectorlEvectorlEfloatgRsPgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEvectorlEfloatgRsPgRsPgR(void *p = 0);
   static void *newArray_vectorlEvectorlEvectorlEfloatgRsPgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEvectorlEfloatgRsPgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEvectorlEfloatgRsPgRsPgR(void *p);
   static void destruct_vectorlEvectorlEvectorlEfloatgRsPgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<vector<float> > >*)
   {
      vector<vector<vector<float> > > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<vector<float> > >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<vector<float> > >", -2, "vector", 210,
                  typeid(vector<vector<vector<float> > >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEvectorlEfloatgRsPgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<vector<float> > >) );
      instance.SetNew(&new_vectorlEvectorlEvectorlEfloatgRsPgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEvectorlEfloatgRsPgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEvectorlEfloatgRsPgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEvectorlEfloatgRsPgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEvectorlEfloatgRsPgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<vector<float> > > >()));

      ::ROOT::AddClassAlternate("vector<vector<vector<float> > >","std::vector<std::vector<std::vector<float, std::allocator<float> >, std::allocator<std::vector<float, std::allocator<float> > > >, std::allocator<std::vector<std::vector<float, std::allocator<float> >, std::allocator<std::vector<float, std::allocator<float> > > > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<vector<float> > >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEvectorlEfloatgRsPgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<vector<float> > >*)0x0)->GetClass();
      vectorlEvectorlEvectorlEfloatgRsPgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEvectorlEfloatgRsPgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEvectorlEfloatgRsPgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<vector<float> > > : new vector<vector<vector<float> > >;
   }
   static void *newArray_vectorlEvectorlEvectorlEfloatgRsPgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<vector<float> > >[nElements] : new vector<vector<vector<float> > >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEvectorlEfloatgRsPgRsPgR(void *p) {
      delete ((vector<vector<vector<float> > >*)p);
   }
   static void deleteArray_vectorlEvectorlEvectorlEfloatgRsPgRsPgR(void *p) {
      delete [] ((vector<vector<vector<float> > >*)p);
   }
   static void destruct_vectorlEvectorlEvectorlEfloatgRsPgRsPgR(void *p) {
      typedef vector<vector<vector<float> > > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<vector<float> > >

namespace ROOT {
   static TClass *vectorlEvectorlEintgRsPgR_Dictionary();
   static void vectorlEvectorlEintgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEintgRsPgR(void *p = 0);
   static void *newArray_vectorlEvectorlEintgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEintgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEintgRsPgR(void *p);
   static void destruct_vectorlEvectorlEintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<int> >*)
   {
      vector<vector<int> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<int> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<int> >", -2, "vector", 210,
                  typeid(vector<vector<int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEintgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<int> >) );
      instance.SetNew(&new_vectorlEvectorlEintgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEintgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEintgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEintgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEintgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<int> > >()));

      ::ROOT::AddClassAlternate("vector<vector<int> >","std::vector<std::vector<int, std::allocator<int> >, std::allocator<std::vector<int, std::allocator<int> > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<int> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<int> >*)0x0)->GetClass();
      vectorlEvectorlEintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEintgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<int> > : new vector<vector<int> >;
   }
   static void *newArray_vectorlEvectorlEintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<int> >[nElements] : new vector<vector<int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEintgRsPgR(void *p) {
      delete ((vector<vector<int> >*)p);
   }
   static void deleteArray_vectorlEvectorlEintgRsPgR(void *p) {
      delete [] ((vector<vector<int> >*)p);
   }
   static void destruct_vectorlEvectorlEintgRsPgR(void *p) {
      typedef vector<vector<int> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<int> >

namespace ROOT {
   static TClass *vectorlEvectorlEfloatgRsPgR_Dictionary();
   static void vectorlEvectorlEfloatgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEfloatgRsPgR(void *p = 0);
   static void *newArray_vectorlEvectorlEfloatgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEfloatgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEfloatgRsPgR(void *p);
   static void destruct_vectorlEvectorlEfloatgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<float> >*)
   {
      vector<vector<float> > *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<float> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<float> >", -2, "vector", 210,
                  typeid(vector<vector<float> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEfloatgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<float> >) );
      instance.SetNew(&new_vectorlEvectorlEfloatgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEfloatgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEfloatgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEfloatgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEfloatgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<float> > >()));

      ::ROOT::AddClassAlternate("vector<vector<float> >","std::vector<std::vector<float, std::allocator<float> >, std::allocator<std::vector<float, std::allocator<float> > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<float> >*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEfloatgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<float> >*)0x0)->GetClass();
      vectorlEvectorlEfloatgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEfloatgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEfloatgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<float> > : new vector<vector<float> >;
   }
   static void *newArray_vectorlEvectorlEfloatgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<float> >[nElements] : new vector<vector<float> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEfloatgRsPgR(void *p) {
      delete ((vector<vector<float> >*)p);
   }
   static void deleteArray_vectorlEvectorlEfloatgRsPgR(void *p) {
      delete [] ((vector<vector<float> >*)p);
   }
   static void destruct_vectorlEvectorlEfloatgRsPgR(void *p) {
      typedef vector<vector<float> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<float> >

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = 0);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 210,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      ::ROOT::AddClassAlternate("vector<int>","std::vector<int, std::allocator<int> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)0x0)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = 0);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 210,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));

      ::ROOT::AddClassAlternate("vector<float>","std::vector<float, std::allocator<float> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<float>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<float>*)0x0)->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete ((vector<float>*)p);
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] ((vector<float>*)p);
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<float>

namespace {
  void TriggerDictionaryInitialization_mydict_Impl() {
    static const char* headers[] = {
"./interface/SetTDRStyle.h",
"./interface/AnalysisUtils.h",
"./interface/FitUtils.h",
"./interface/Na22SpectrumAnalyzer.h",
0
    };
    static const char* includePaths[] = {
"/home/cmsdaq/alio_guglielmi_corr/Lab5015Analysis",
"/usr/include/root",
"/usr/include/root",
"/home/cmsdaq/alio_guglielmi_corr/Lab5015Analysis/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "mydict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
class __attribute__((annotate("$clingAutoload$./interface/AnalysisUtils.h")))  EventClass;
class __attribute__((annotate("$clingAutoload$./interface/AnalysisUtils.h")))  ModuleEventClass;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "mydict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "./interface/SetTDRStyle.h"
#include "./interface/AnalysisUtils.h"
#include "./interface/FitUtils.h"
#include "./interface/Na22SpectrumAnalyzer.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"EventClass", payloadCode, "@",
"ModuleEventClass", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("mydict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_mydict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_mydict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_mydict() {
  TriggerDictionaryInitialization_mydict_Impl();
}
