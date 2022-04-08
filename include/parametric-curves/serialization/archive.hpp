#ifndef __parametric_curves_serialization_archive_hpp__
#define __parametric_curves_serialization_archive_hpp__

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/xml_iarchive.hpp>
#include <boost/archive/xml_oarchive.hpp>
#include <fstream>
#include <stdexcept>
#include <string>

namespace serialization {
template <class Derived>
struct Serializable {
 private:
  Derived& derived() { return *static_cast<Derived*>(this); }
  const Derived& derived() const { return *static_cast<const Derived*>(this); }

 public:
  /// \brief Loads a Derived object from a text file.
  void loadFromText(const std::string& filename) {
    std::ifstream ifs(filename.c_str());
    if (ifs) {
      boost::archive::text_iarchive ia(ifs);
      ia >> derived();
    } else {
      const std::string exception_message(filename +
                                          " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }

  /// \brief Saved a Derived object as a text file.
  void saveAsText(const std::string& filename) const {
    std::ofstream ofs(filename.c_str());
    if (ofs) {
      boost::archive::text_oarchive oa(ofs);
      oa << derived();
    } else {
      const std::string exception_message(filename +
                                          " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }

  /// \brief Loads a Derived object from an XML file.
  void loadFromXML(const std::string& filename, const std::string& tag_name) {
    assert(!tag_name.empty());
    std::ifstream ifs(filename.c_str());
    if (ifs) {
      boost::archive::xml_iarchive ia(ifs);
      ia >> boost::serialization::make_nvp(tag_name.c_str(), derived());
    } else {
      const std::string exception_message(filename +
                                          " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }

  /// \brief Saved a Derived object as an XML file.
  void saveAsXML(const std::string& filename,
                 const std::string& tag_name) const {
    assert(!tag_name.empty());
    std::ofstream ofs(filename.c_str());
    if (ofs) {
      boost::archive::xml_oarchive oa(ofs);
      oa << boost::serialization::make_nvp(tag_name.c_str(), derived());
    } else {
      const std::string exception_message(filename +
                                          " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }

  /// \brief Loads a Derived object from an binary file.
  void loadFromBinary(const std::string& filename) {
    std::ifstream ifs(filename.c_str());
    if (ifs) {
      boost::archive::binary_iarchive ia(ifs);
      ia >> derived();
    } else {
      const std::string exception_message(filename +
                                          " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }

  /// \brief Saved a Derived object as an binary file.
  void saveAsBinary(const std::string& filename) const {
    std::ofstream ofs(filename.c_str());
    if (ofs) {
      boost::archive::binary_oarchive oa(ofs);
      oa << derived();
    } else {
      const std::string exception_message(filename +
                                          " does not seem to be a valid file.");
      throw std::invalid_argument(exception_message);
    }
  }
};

}  // namespace serialization

#endif  // ifndef __parametric_curves_serialization_archive_hpp__
