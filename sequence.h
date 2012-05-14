#ifndef SEQUENCE_H
#define SEQUENCE_H

#include <memory>

//----------------------------------------------------------------

template <typename T>
class sequence {
public:
	typedef std::shared_ptr<sequence<T> > ptr;

	virtual ~sequence() {}
	virtual T operator()() = 0;
};

//----------------------------------------------------------------

#endif
