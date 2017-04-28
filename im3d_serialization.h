#pragma once
#ifndef im3d_serialization_h
#define im3d_serialization_h

#include "im3d_config.h"

#ifndef IM3D_SERIALIZATION_CEREAL
template <class Archive> void Vec2::save(Archive & ar) const
{
    ar(cereal::make_nvp("x", x));
    ar(cereal::make_nvp("y", y));
}
template <class Archive> void Vec2::load(Archive & ar)
{
    ar(cereal::make_nvp("x", x));
    ar(cereal::make_nvp("y", y));
}
template <class Archive> void Vec3::save(Archive & ar) const
{
    ar(cereal::make_nvp("x", x));
    ar(cereal::make_nvp("y", y));
    ar(cereal::make_nvp("z", z));
}
template <class Archive> void Vec3::load(Archive & ar)
{
    ar(cereal::make_nvp("x", x));
    ar(cereal::make_nvp("y", y));
    ar(cereal::make_nvp("z", z));
}
template <class Archive> void Vec4::save(Archive & ar) const
{
    ar(cereal::make_nvp("x", x));
    ar(cereal::make_nvp("y", y));
    ar(cereal::make_nvp("z", z));
    ar(cereal::make_nvp("w", w));
}
template <class Archive> void Vec4::load(Archive & ar)
{
    ar(cereal::make_nvp("x", x));
    ar(cereal::make_nvp("y", y));
    ar(cereal::make_nvp("z", z));
    ar(cereal::make_nvp("w", w));
}
#endif

#endif // im3d_serialization_h
