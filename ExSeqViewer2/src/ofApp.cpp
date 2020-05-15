#include "ofApp.h"

//For loading the CSV, set a few global paramters
#define CSV_COL_X 0
#define CSV_COL_Y 1
#define CSV_COL_Z 2
#define CSV_COL_GENESYMBOL 3
#define CSV_COL_READTYPE 4
#define CSV_COL_DIDALIGN 5
#define CSV_COL_INMORPH 6
//#define CSV_COL_CELLNUMBER 6

//For handling different readtypes, create a color mapping
#define PTYPE_INTRON 0
#define PTYPE_EXON 1
#define PTYPE_INTERGENIC 2
#define PTYPE_RRNA18s 3
#define PTYPE_RRNA28s 4
#define PTYPE_RRNA45s 5
#define PTYPE_NOALIGNMENT 6
#define PTYPE_NCRNA 7
#define PTYPE_TFMRNA 8
#define PTYPE_TRNA 9

//--------------------------------------------------------------
void ofApp::setup(){
    ofSetVerticalSync(true);
    //      ofEnableSmoothing(); //Not sure that this does

    verdana.load("verdana.ttf", 15, true, true);
    verdana.setLineHeight(34.0f);
    verdana.setLetterSpacing(1.035);
    //Setup gui
    gui.setup(); // most of the time you don't need a name but don't forget to call setup
    gui.setPosition(28, 150);
    gui.add(inMorphologyOnly.setup("Reads in annotated cells only", true));
    gui.add(showAllLabels.setup("Show all labels", true));
    gui.add(cursorTrack.setup("Show labels under mouse", true));
    gui.add(singleNeuronOnly.setup("Show only one neuron", false));
    gui.add(singleNeuronIndex.setup("neuron index", 3, 1, 13));
    gui.add(alphaNeurons.setup("alphaNeurons", 40, 0, 255));
    
    gui.add(punctaSize.setup("punctaSize", 10, 0, 50));
    gui.add(fontSize.setup("fontSize", 10, 0, 50));
    showGUI = true;
    
    //Initilize a set of colors
    colorMapForPunctaType.clear();
    colorMapForPunctaType.push_back(ofColor(213,66.,213.,250.));     //intron magenta
    colorMapForPunctaType.push_back(ofColor(66.,6.,213., 250.));   //exon blue
    colorMapForPunctaType.push_back(ofColor(0.,0,0,250.));          //intergenic cyan shouldn't be any of these
    colorMapForPunctaType.push_back(ofColor(255,255,255,230.));     //rrna-18s white
    colorMapForPunctaType.push_back(ofColor(50.,125,125,250.));  //rrna-28s cyan
    colorMapForPunctaType.push_back(ofColor(0.,0.,0.,230.));        //rrna-45s black
    colorMapForPunctaType.push_back(ofColor(100.,100,100,230.));    //noalignment gray
    colorMapForPunctaType.push_back(ofColor(197,77.,72.,250.));      //ncRNA
    colorMapForPunctaType.push_back(ofColor(213.,140.,66.,250.));   //transcription factor
    colorMapForPunctaType.push_back(ofColor(233.,166.,30.,250.));   //tRNA
    
    //Load all the neuron meshes, if available
    string filename;
    for(int fov_idx = 8; fov_idx < 11; fov_idx += 1) {
        for(int idx = 0; idx < 20; idx += 1) {
            
            filename = "fov_" +std::to_string(fov_idx	) + "_neuron_" + std::to_string(idx+1)+".ply";
            ofFile file(filename);
            
            if (file.exists()) {
                ofMesh mesh;
                mesh.load(filename);
                this->neuronMeshes.push_back(mesh);
            }
        }
    }

    //Load all the nucleus meshes, if available
//    for(int fov_idx = 8; fov_idx < 11; fov_idx += 1) {
//        for(int idx = 0; idx < 20; idx += 1) {
//
//            filename = "fov_" +std::to_string(fov_idx    ) + "_nucleus_" + std::to_string(idx+1)+".ply";
//            ofFile file(filename);
//
//            if (file.exists()) {
//                ofMesh mesh;
//                mesh.load(filename);
//                this->nucleusMeshes.push_back(mesh);
//            }
//        }
//    }

    
    
    meshpuncta.setMode(OF_PRIMITIVE_POINTS);
    ofEnableDepthTest();
    
    ofSetVerticalSync(true);
    
    
    // Clear the punctalist just to be safe
    punctalist.clear();
    // Loop through mulitple fovs worth of puncta
    for(int fov_idx = 8; fov_idx < 11; fov_idx += 1) {

        filename = "fov_" +std::to_string(fov_idx) + "_puncta.csv";
        ofFile file(filename);
        if (file.exists()) {
            loadCSVIntoPointCloud(filename);
            cout << "Loaded " << filename << ", total puncta =" <<  punctalist.size() << endl;
        }
    }
   
    
    makePunctaMesh();

//    image.loadImage("Slice_4x_final.bmp");
//    image.load("Slice_4x_final.bmp");
//    image.resize
    
    
    ofVec3f pos;
    int N = (int)punctalist.size();
    for(int idx = 0; idx < N; idx += 1) {
        pos = pos+punctalist[idx].pos;
    }
    pos.x = pos.x/N; pos.y=pos.y/N; pos.z= pos.z/N;
    
    camera.setPosition(pos.x, pos.y, pos.z + 400);
    camera.lookAt(pos);
    camera.rotate(270, ofVec3f(0,0,1));
//    plane.setPosition(pos.y/N, pos.x/N, pos.z/N);
    
    cout << pos << endl;
    cout << camera.getGlobalPosition() << endl;
}

//--------------------------------------------------------------
void ofApp::update(){
 
    // Check if the toggled has been thrown
    if ( last_inMorphologyOnly_value!=inMorphologyOnly ||
        last_singleNeuronIndex_value != singleNeuronIndex ||
        last_singleNeuronOnly_value != singleNeuronOnly){
        
        makePunctaMesh();
        
        last_inMorphologyOnly_value = inMorphologyOnly;
        last_singleNeuronIndex_value = singleNeuronIndex;
        last_singleNeuronOnly_value = singleNeuronOnly;
    }
    
    
}

//--------------------------------------------------------------
void ofApp::draw(){
    camera.setFarClip(10000.);
    
    ofBackgroundGradient(ofColor(255), ofColor(255));
    
//    image.getTextureReference().bind();
//    plane.draw();
//    image.getTextureReference().unbind();
    
    
    ofSetColor(255);
    

    camera.begin();
    
    

    ofSetColor(ofColor(255.,255.,0.), 120  );
    glPointSize(punctaSize); // default 10 make the points bigger
    meshpuncta.drawVertices();
    
    
//    ofSetColor(ofColor(0,255.,0.), 20  );
//    mesh.drawWireframe();
//    mesh.drawFaces();
//    glPointSize(2);

    // Setting the color for the morphology
    ofSetColor(ofColor(0,200.,0.), alphaNeurons  );
    for(int idx = 0; idx < neuronMeshes.size(); idx += 1) {
        
        //If we're set to only show one cell at a time
        //note that there is a +1 for the difference in 0-indexed vs 1-indexing
        if (singleNeuronOnly && singleNeuronIndex!=(idx+1)){
            continue;
        }
        
        neuronMeshes[idx].drawWireframe();
    }
    
//    ofSetColor(ofColor(0,0,200.), 255  );
//    for(int idx = 0; idx < nucleusMeshes.size(); idx += 1) {
//        nucleusMeshes[idx].drawWireframe();
//    }
    
    glPointSize(2);
    

    
//    ofSetColor(ofColor::white);
    
    
    
    camera.end();
    
    int n = (int) meshpuncta.getNumVertices();
    float nearestDistance = 0;
    ofVec2f nearestVertex;
    int nearestIndex = 0;
    ofVec2f mouse(mouseX, mouseY);
    ofVec2f pos2labelshow;
    for(int i = 0; i < n; i++) {

        ofVec3f cur = camera.worldToScreen(meshpuncta.getVertex(i));
        float distance = cur.distance(mouse);
        if(i == 0 || distance < nearestDistance) {
            
            nearestDistance = distance;
            nearestVertex = cur;
            nearestIndex = i;
        }

        if(punctalist[punctaShown[i]].readtype== PTYPE_EXON ||
           punctalist[punctaShown[i]].readtype== PTYPE_INTRON ||
           punctalist[punctaShown[i]].readtype== PTYPE_NCRNA ||
           punctalist[punctaShown[i]].readtype== PTYPE_TFMRNA
            )
        {
            ofNoFill();
            ofSetColor(punctalist[punctaShown[i]].color);
            ofSetLineWidth(2);
            ofDrawCircle(cur, punctaSize/2+1);
            ofSetLineWidth(1);
            
            if(showAllLabels){
                ofVec2f offset(2*punctaSize, -1*punctaSize);
                pos2labelshow = cur + offset;

                //To show the gene neames in the old fasion
                //                ofDrawBitmapStringHighlight(punctalist[punctaShown[i]].gene_symbol, cur + offset);
                verdana.drawString(punctalist[punctaShown[i]].gene_symbol, pos2labelshow.x, pos2labelshow.y);
                
                ofSetLineWidth(2);
                ofDrawLine(cur, pos2labelshow);
                
                
            }
        }
        else //if (punctalist[punctaShown[i]].readtype== PTYPE_RRNA18s ||
             //   punctalist[punctaShown[i]].readtype== PTYPE_RRNA28s ||
             //   punctalist[punctaShown[i]].readtype== PTYPE_RRNA45s)
        {
            ofNoFill();
            ofSetColor(ofColor::black);
            ofSetLineWidth(1);
            ofDrawCircle(cur, punctaSize/2+1);
            
        }
        
    }
    
    if (cursorTrack){
        ofSetColor(ofColor::gray);
        ofDrawLine(nearestVertex, mouse);
        
        ofNoFill();
        ofSetColor(ofColor::black);
        ofSetLineWidth(2);
        ofDrawCircle(nearestVertex, punctaSize/2+1);
        ofSetLineWidth(1);
        
        ofVec2f offset(punctaSize, -1*punctaSize);
    //    ofDrawBitmapStringHighlight(ofToString(nearestIndex), mouse + offset);
        //Punctashown is an additional vector to keep track of the subset of punctalist that is current
        //Being shown. It is indexed just as the the mesh points are, which can be a subset of punctalist.
        ofDrawBitmapStringHighlight(punctalist[punctaShown[nearestIndex]].gene_symbol, mouse + offset);
    }
    
    
    //--
    //draw 2d text notes on top
    //--
    //
    if(showGUI){
        ofPushStyle();
        {
            ofSetColor(0, 0, 0);
            
            ofDrawBitmapStringHighlight("ExSeqViewer Beta", 30, 30, ofColor::darkSlateGray, ofColor::white);
            
            stringstream message;
            message << "Navigation keys:" << endl;
            message << "Click and drag to orbit" << endl;
            message << "Control+click drag with right mouse to zoom" << endl;
            message << "Press 'c' to specify a rotation origin" << endl;
            message << "Hold 'h' and drag with left mouse to pan.'x' hides menu" << endl;
            message << "" << endl;
            ofDrawBitmapStringHighlight(message.str(), 30, 60, ofColor::darkSlateGray, ofColor::white);
            
            ofDrawBitmapStringHighlight("Exon", 30,130,ofColor::darkSlateGray, colorMapForPunctaType[PTYPE_EXON]);
            ofDrawBitmapStringHighlight("Intron", 90,130,ofColor::darkSlateGray, colorMapForPunctaType[PTYPE_INTRON]);
//            ofDrawBitmapStringHighlight("rRNA", 160,130,ofColor::darkSlateGray, colorMapForPunctaType[PTYPE_RRNA]);
            
        }
        ofPopStyle();
        //Render the GUI last so it's on top of the 3D render
        gui.draw();
        
    }
    

}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
 
    if (key == 'c') {
        camera.toggleCursorDrawEnabled();
    }
    if (key == 'u') {
        camera.toggleFixUpDirectionEnabled();
    }
    if (key == 'x') {
        if(showGUI)
            showGUI = false;
        else
            showGUI = true;
    }
    // try reloading the font to see if we can rescale
    if (key == 'f') {
        verdana.load("verdana.ttf", fontSize, true, true);
        verdana.setLineHeight(34.0f);
        verdana.setLetterSpacing(1.035);
    }
    
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){
    
}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y){
    
}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){
    
}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){
    
}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){
    
}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){
    
}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){
    
}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){
    
}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){
    
}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 
    
}
//--------------------------------------------------------------

void ofApp::loadCSVIntoPointCloud(const string &path){
    
   
    
    csv.load(path);
    

    float x = 0.0;
    float y = 0.0;
    float z = 0.0;
    
    float r = 0.0;
    float g = 0.0;
    float b = 0.0;
    float a = 0.0;
    
    float val = 0.;
    
    int puncta_list_start_idx = (int)punctalist.size();
    int insert_idx;
    for(int idx = 0; idx < csv.getNumRows(); idx += 1) {
        
        //Create a new puncta object and push it to the back
        punctalist.push_back(Puncta());
        
        x = ::atof(csv[idx][CSV_COL_X].c_str());
        y = ::atof(csv[idx][CSV_COL_Y].c_str());
        z = ::atof(csv[idx][CSV_COL_Z].c_str());
        
        ofVec3f pos(x,y,z);
        
        insert_idx = (int)punctalist.size()-1;
        
        punctalist[insert_idx].pos = pos;
        
        punctalist[insert_idx].readtype = :: atof(csv[idx][CSV_COL_READTYPE].c_str())-1;
        punctalist[insert_idx].color = colorMapForPunctaType[punctalist[insert_idx].readtype];
        
        punctalist[insert_idx].gene_symbol = csv[idx][CSV_COL_GENESYMBOL];
        punctalist[insert_idx].doShow = true;
        punctalist[insert_idx].aligned= :: atof(csv[idx][CSV_COL_DIDALIGN].c_str());
        
        punctalist[insert_idx].inMorphology = :: atof(csv[idx][CSV_COL_INMORPH].c_str())>0;
        punctalist[insert_idx].cellnumber = :: atof(csv[idx][CSV_COL_INMORPH].c_str());
        
        
    }
    
    

    
    
    cout << "Done..\n";
}

//--------------------------------------------------------------
void ofApp::makePunctaMesh()
{
    meshpuncta.clear();
    punctaShown.clear();
    for(int idx = 0; idx < punctalist.size(); idx += 1) {
        
        // create a vector punctashown, that keeps track of the puncta
        if (inMorphologyOnly && !punctalist[idx].inMorphology )
            continue;
        
        if (singleNeuronOnly && singleNeuronIndex!=punctalist[idx].cellnumber){
            continue;
        }
        
        meshpuncta.addColor(punctalist[idx].color);
        meshpuncta.addVertex(punctalist[idx].pos);
        punctaShown.push_back(idx);
    }
}
