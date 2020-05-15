//
//  ofPuncta.hpp
//  ExSeqViewer2
//
//  Created by Daniel Goodwin on 11/2/18.
//

#ifndef ofPuncta_hpp
#define ofPuncta_hpp

#include <stdio.h>
#include "ofMain.h"

class ofPuncta {
//private:
    // vars
    ofVec3f pos;
    string gene_symbol;
    string readtype;
    bool aligned;
    ofColor color;
//    float x;
//    float y;
//    float z;
//    float speedX;
//    float speedY;
    
//public:
    // methods
    void update();
    void draw();
    
    // constructor
    ofPuncta(float _x, float _y, float _z);
};

#endif /* ofPuncta_hpp */
