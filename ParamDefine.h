#pragma once
#ifndef PARAMDEFINE_H
#define PARAMDEFINE_H

/////////////////////////////////////////////////
// STL //////////////////////////////////////////
#include <iostream>
#include <memory>

/////////////////////////////////////////////////
// VTK //////////////////////////////////////////
#include "vtkSmartPointer.h"

#define STEP_SEPARATOR 0xffffffff // the largest unsigned int in 32bit machine
#define CONTACT_SEPARATOR 0xfffffffe

#define NODE_TRIANGLE 0x0001
#define TRIANGLE_NODE 0x0002
#define EDGE_EDGE	  0x0004
#define NODE_EDGE	  0x0008
#define EDGE_NODE	  0x0010
#define NODE_NODE	  0x0020

#define DEFAULT_TIMERCALLBACK TimerCallback
#define DEFAULT_KEYPRESSCALLBACK KeypressCallback
#define DEFAULT_WINDOWCALLBACK WindowModifiedCallback

const unsigned char red[3] = { 255, 0, 0 };
const unsigned char green[3] = { 0, 255, 0 };
const unsigned char blue[3] = { 0, 0, 255 };
const unsigned char cyan[3] = { 0, 255, 255 };

#endif