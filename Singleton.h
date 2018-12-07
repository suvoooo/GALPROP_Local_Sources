#ifndef _utl_Singleton_h_
#define _utl_Singleton_h_

namespace utl {

  /**
   * \class Singleton Singleton.h utl/Singleton.h
   *
   * \brief Curiously Recurring Template Pattern (CRTP) for Meyers singleton
   *
   * The singleton class is implemented as follows
   * \code
   * #include <utl/Singleton.h>
   *
   * class SomeClass : public utl::Singleton<SomeClass> {
   *   ...
   * private:
   *   // prevent creation, destruction
   *   SomeClass() { }
   *   ~SomeClass() { }
   *
   *   friend class utl::Singleton<SomeClass>;
   * };
   * \endcode
   * Singleton automatically prevents copying of the derived class.
   *
   * \author Darko Veberic
   * \date 9 Aug 2006
   * \version $Id: Singleton.h 7249 2008-05-07 22:03:06Z darko $
   * \ingroup stl
   */

  template<typename T>
  class Singleton {
  public:
    static
    T&
    GetInstance()
#ifdef __MAKECINT__
    ;
#else
    {
      static T instance;
      return instance;
    }
#endif

  protected:
    // derived class can call ctor and dtor
    Singleton() { }
    ~Singleton() { }

  private:
    // no one should do copies
    Singleton(const Singleton&);
    Singleton& operator=(const Singleton&);

  };


  /**
   * \class LeakingSingleton Singleton.h utl/Singleton.h
   *
   * \brief CRTP for leaking singleton
   *
   * This type of creation (Gamma singleton) leaks the object at
   * the end of the run, i.e. class destructor does not get called
   * in at_exit().
   *
   * This singleton can be implemented as follows
   * \code
   * #include <utl/Singleton.h>
   *
   * class SomeClass : public utl::LeakingSingleton<SomeClass> {
   *   ...
   * private:
   *   // prevent creation, destruction
   *   SomeClass() { }
   *   ~SomeClass() { }
   *
   *   friend class utl::LeakingSingleton<SomeClass>;
   * };
   * \endcode
   * LeakingSingleton automatically prevents copying of the derived
   * class.
   *
   * \author Darko Veberic
   * \date 9 Aug 2006
   * \version $Id: Singleton.h 7249 2008-05-07 22:03:06Z darko $
   * \ingroup stl
   */

  template<class T>
  class LeakingSingleton {
  public:
    static
    T&
    GetInstance()
    {
      if (!fgInstance)
        fgInstance = new T;
      return *fgInstance;
    }

  protected:
    // derived class can call ctor and dtor
    LeakingSingleton() { }
    // will never get called anyway
    ~LeakingSingleton() { }

  private:
    // no one should do copies
    LeakingSingleton(const LeakingSingleton&);
    LeakingSingleton& operator=(const LeakingSingleton&);

    static T* fgInstance;
  };

  template<class T>
  T* LeakingSingleton<T>::fgInstance = 0;

}


#endif
