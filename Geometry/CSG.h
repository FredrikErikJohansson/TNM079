
/*************************************************************************************************
 *
 * Modeling and animation (TNM079) 2007
 * Code base for lab assignments. Copyright:
 *   Gunnar Johansson (gunnar.johansson@itn.liu.se)
 *   Ken Museth (ken.museth@itn.liu.se)
 *   Michael Bang Nielsen (bang@daimi.au.dk)
 *   Ola Nilsson (ola.nilsson@itn.liu.se)
 *   Andreas Sderstrm (andreas.soderstrom@itn.liu.se)
 *
 *************************************************************************************************/
#ifndef __CSG_H__
#define __CSG_H__

#include "Geometry/Implicit.h"
#include <cmath>

/*! \brief CSG Operator base class */
class CSG_Operator : public Implicit {

protected:
  //! Constructor
  CSG_Operator(Implicit *l, Implicit *r) : left(l), right(r) {}

  //! Pointers to left and right child nodes
  Implicit *left, *right;
};

/*! \brief Union boolean operation */
class Union : public CSG_Operator {
public:
  Union(Implicit *l, Implicit *r) : CSG_Operator(l, r) {
    // Compute the resulting (axis aligned) bounding box from
    // the left and right children
    mBox = BoxUnion(l->GetBoundingBox(), r->GetBoundingBox());
  }

  virtual float GetValue(float x, float y, float z) const {
    TransformW2O(x, y, z);
    return (std::min(left->GetValue(x,y,z), right->GetValue(x,y,z)));
  }
};

/*! \brief Intersection boolean operation */
class Intersection : public CSG_Operator {
public:
  Intersection(Implicit *l, Implicit *r) : CSG_Operator(l, r) {
    mBox = BoxIntersection(l->GetBoundingBox(), r->GetBoundingBox());
  }

  virtual float GetValue(float x, float y, float z) const { 
    TransformW2O(x, y, z);
    return (std::max(left->GetValue(x,y,z), right->GetValue(x,y,z)));
  }
};

/*! \brief Difference boolean operation */
class Difference : public CSG_Operator {
public:
  Difference(Implicit *l, Implicit *r) : CSG_Operator(l, r) {
    mBox = l->GetBoundingBox();
  }

  virtual float GetValue(float x, float y, float z) const {
    TransformW2O(x, y, z);
    return (std::max(left->GetValue(x,y,z), -right->GetValue(x,y,z)));
  }
};

/*! \brief BlendedUnion boolean operation */
class BlendedUnion : public CSG_Operator {
public:
  BlendedUnion(Implicit *l, Implicit *r, int blend)
      : CSG_Operator(l, r), mBlend(blend) {
    mBox = BoxUnion(l->GetBoundingBox(), r->GetBoundingBox());
  }

  virtual float GetValue(float x, float y, float z) const {
    TransformW2O(x, y, z);
    float Da = exp(-left->GetValue(x, y, z));
    float Db = exp(-right->GetValue(x, y, z));
    float val = std::pow((std::pow(Da, mBlend) + std::pow(Db, mBlend)), (1.0f / mBlend));
  	return 1 - val;
  }

protected:
  int mBlend;
};

/*! \brief BlendedIntersection boolean operation */
class BlendedIntersection : public CSG_Operator {
public:
  BlendedIntersection(Implicit *l, Implicit *r, int blend)
      : CSG_Operator(l, r), mBlend(blend) {
    mBox = BoxUnion(l->GetBoundingBox(), r->GetBoundingBox());
  }

  virtual float GetValue(float x, float y, float z) const {
      TransformW2O(x, y, z);
      float Da = exp(-left->GetValue(x, y, z));
      float Db = exp(-right->GetValue(x, y, z));
      float val = std::pow((std::pow(Da, -mBlend) + std::pow(Db, -mBlend)), (-1.0f / mBlend));
      return 1 - val;
  }

protected:
  int mBlend;
};

/*! \brief BlendedDifference boolean operation */
class BlendedDifference : public CSG_Operator {
public:
  BlendedDifference(Implicit *l, Implicit *r, int blend)
      : CSG_Operator(l, r), mBlend(blend) {
    mBox = l->GetBoundingBox();
  }

  virtual float GetValue(float x, float y, float z) const {
      TransformW2O(x, y, z);
      float Da = exp(-left->GetValue(x, y, z));
      float Db = exp(right->GetValue(x, y, z));
      float val = std::pow((std::pow(Da, mBlend) + std::pow(Db, mBlend)), (1.0f / mBlend));
      return 1 - val;
  }

protected:
  int mBlend;
};

#endif
