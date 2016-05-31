//
//  CircularBuffer.hpp
//  mtg_new
//
//  Created by Guillaume Rizk 
//  Copyright © 2016 Guillaume Rizk. All rights reserved.
//

#ifndef CircularBuffer_hpp
#define CircularBuffer_hpp

#include <stdio.h>


#include <stdlib.h>



template <typename elem_t>  class CircularBuffer
{

	friend class itCB;
	
public:
	
	//constructor, len in power of two  : 2^powlen
	CircularBuffer(int powlen) : _idx(this)
	{
		_tai   = (1LL << powlen);
		_buffer = (elem_t* ) calloc(_tai,sizeof(elem_t));
		_mask =  _tai -1;
	}
	
	//def construc
	CircularBuffer() : _idx(this)
	{
		int powlen = 10;
		_tai   = (1LL << powlen);
		_buffer = (elem_t* ) calloc(_tai,sizeof(elem_t));
		_mask =  _tai -1;
	}
	
	void resize(int powlen)
	{
		_tai   = (1LL << powlen);
		_buffer = (elem_t* ) realloc(_buffer,sizeof(elem_t)*_tai);
		_mask =  _tai -1;
	}
	
	
	void push(elem_t new_elem)
	{
		_idx() = new_elem ;
		//_buffer[_idx] = new_elem;
		_idx++;
	}
	
	~CircularBuffer()
	{
		free(_buffer);
	};
	
	void clear()
	{
		memset(_buffer,0,sizeof(elem_t)*_tai);
		_idx = itCB(this);
	}
	

	class itCB // iterator of circularbuffer
	{
	public:
		u_int64_t _idx;
		CircularBuffer<elem_t> * _ref;
		
		itCB (CircularBuffer *ref)
		{
			_idx = 0;
			_ref=ref;
		}
		
		elem_t & operator()()
		{
			return _ref->_buffer[_idx];
		}
		elem_t & item()
		{
			return _ref->_buffer[_idx];
		}
		
		void set(u_int64_t val)
		{
			_idx = val & _ref->_mask;
		}
		
		itCB   operator+( u_int64_t rhs) //member func
		{
			itCB  nit =  itCB(this->_ref);
			nit._idx =  (this->_idx + rhs) & _ref->_mask;
			return nit;
		}
		itCB   operator-( u_int64_t rhs) //member func
		{
			itCB  nit =  itCB(this->_ref);
			nit._idx =  (this->_idx - rhs) & _ref->_mask;
			return nit;
		}
		void operator++(int) //postfix operator
		{
			_idx = (_idx+1) & _ref->_mask;
		}
		void operator--(int)
		{
			_idx = (_idx-1) & _ref->_mask;
		}
	};
	
private:
	
	elem_t * _buffer;
	u_int64_t _tai;
	u_int64_t _mask;
	itCB _idx; //index to next free cell
	
};




#endif /* CircularBuffer_hpp */

