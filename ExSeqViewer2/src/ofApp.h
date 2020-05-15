#pragma once

#include "ofMain.h"
#include "ofxCsv.h"
#include "ofxGui.h"
#include "ofxGrabCam.h"

class ofApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);
		
    
    //Variables that will be useful
//        ofMesh mesh;
        ofMesh meshpuncta;
//        ofEasyCam cam;
        ofxGrabCam camera;
        ofxCsv csv;
    
    
    bool showGUI;
    
    
    //A class to hold the puncta information in memory
    struct Puncta {
        ofVec3f pos;
        string gene_symbol;
        int readtype;
        bool aligned;
        ofColor color;
        bool doShow;
        bool inMorphology;
        int cellnumber;
    };
    //Collection of all the puncta with their associated metadata
    vector<Puncta> punctalist;
    //Create a handcoded set of colors for the puncta by readtypes
    vector<ofColor> colorMapForPunctaType;
    //Create a set of meshes for each piece of morphology
    vector<ofMesh> neuronMeshes;
    
//    vector<ofMesh> nucleusMeshes;
    
    void loadCSVIntoPointCloud(const string &path);
    
    void makePunctaMesh();
    
    
    //GUI variables
    ofxPanel gui;
    ofxToggle inMorphologyOnly;
    ofxToggle showAllLabels;
    ofxToggle cursorTrack;
    ofxToggle singleNeuronOnly;
    
    bool last_inMorphologyOnly_value;
    bool last_singleNeuronOnly_value;
    int last_singleNeuronIndex_value;
    
    vector<int> punctaShown;
    
    ofVec3f data_centroid;
    ofxFloatSlider alphaNeurons;
    ofxFloatSlider punctaSize;
    ofxFloatSlider fontSize;
    ofxIntSlider singleNeuronIndex;
    
    
    
//    void alphaNeuronsChanged(float & alphaNeuronsChanged);

//    ofPlanePrimitive plane;
    ofTrueTypeFont    verdana;
    
    ofImage image;
    
};
