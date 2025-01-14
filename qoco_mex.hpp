#ifndef QOCO_MEX_HPP
#define QOCO_MEX_HPP
#include <stdint.h>
#include <string>
#include <cstring>
#include <typeinfo>
#include <math.h>
#include "qoco.h"
#include "mex.h"

#define QOCO_MEX_SIGNATURE 0x271C1A7A
template <class base>
class qoco_mex_handle
{
public:
    qoco_mex_handle(base *ptr) : m_ptr(ptr), m_name(typeid(base).name()) { m_signature = QOCO_MEX_SIGNATURE; }
    ~qoco_mex_handle()
    {
        m_signature = 0;
        delete m_ptr;
    }
    bool isValid() { return ((m_signature == QOCO_MEX_SIGNATURE) && !strcmp(m_name.c_str(), typeid(base).name())); }
    base *ptr() { return m_ptr; }

private:
    uint32_t m_signature;
    std::string m_name;
    base *m_ptr;
};

template <class base>
inline mxArray *convertPtr2Mat(base *ptr)
{
    mexLock();
    mxArray *out = mxCreateNumericMatrix(1, 1, mxUINT64_CLASS, mxREAL);
    *((uint64_t *)mxGetData(out)) = reinterpret_cast<uint64_t>(new qoco_mex_handle<base>(ptr));
    return out;
}

template <class base>
inline qoco_mex_handle<base> *convertMat2HandlePtr(const mxArray *in)
{
    if (mxGetNumberOfElements(in) != 1 || mxGetClassID(in) != mxUINT64_CLASS || mxIsComplex(in))
        mexErrMsgTxt("Input must be a real uint64 scalar.");
    qoco_mex_handle<base> *ptr = reinterpret_cast<qoco_mex_handle<base> *>(*((uint64_t *)mxGetData(in)));
    if (!ptr->isValid())
        mexErrMsgTxt("Handle not valid.");
    return ptr;
}

template <class base>
inline base *convertMat2Ptr(const mxArray *in)
{
    return convertMat2HandlePtr<base>(in)->ptr();
}

template <class base>
inline void destroyObject(const mxArray *in)
{
    delete convertMat2HandlePtr<base>(in);
    mexUnlock();
}

#endif // QOCO_MEX_HPP