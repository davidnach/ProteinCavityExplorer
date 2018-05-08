
declare var pv:any;

import {UserInputService} from '../userinput.service';
import { Component, OnInit, } from '@angular/core';
import { GlobalData} from '../globaldata.service';

@Component({
  selector: 'app-experiment-results',
  templateUrl: './experiment-results.component.html',
  styleUrls: ['./experiment-results.component.css']
})
export class ExperimentResultsComponent   {
    pid : string;
    experimentId : string;
    color;
    viewer;
    parent;
    pocketNums;
    pocketNumber:number;


  constructor(private userInputService : UserInputService,
private globaldata : GlobalData) {
	//this.pocketNums = Array.from(this.globaldata.getMap()).keys());
  this.pocketNums = [1,2,3,4];

  }


  ngOnInit() {
	  this.experimentId = this.globaldata.getExperimentId();
    this.pid = this.globaldata.getPid();
    console.log(this.pid)
	  this.parent = document.getElementById('viewer');
	  this.color = pv.color;
	  console.log(this.parent);
  }
 ngAfterContentInit() {
	var options = {
          width: 600,
          height: 600,
          antialias: true,
          quality : 'high'
         };

        console.log(this.parent);
 	this.viewer = pv.Viewer(this.parent, options);
	var viewer_ = this.viewer;
 	  this.color = pv.color;
	  console.log(this.viewer);
          var url = "http://www.rcsb.org/pdb/files/";
          url = url + this.pid + ".pdb";
          pv.io.fetchPdb(url, function(structure) {
           viewer_.cartoon('protein', structure);
            viewer_.autoZoom();
            viewer_.centerOn(structure);

	 });
}

tofetchindex(){
  console.log(this.pocketNumber);
}

colorPocketResidues(pocketNum : number){
console.log(pocketNum);
var color_ = this.color;
this.viewer.forEach(function(go) {
    go.colorBy(color_.uniform('red'), go.select({rindices : [1,2,3,4,5,6]}));
     });
    this.viewer.requestRedraw();
}
/*
function setColorForAtom(go, atoms, color) {
    var view = go.structure().createEmptyView();
    var i;
    for(i = 0; i < atoms.length; i++){
	view.addAtom(atoms[i]);
    }
    go.colorBy(pv.color.uniform(color), view);
}

// variable to store the previously picked atom. Required for resetting the color
// whenever the mouse moves.
var prevPicked = null;
// add mouse move event listener to the div element containing the viewer. Whenever
// the mouse moves, use viewer.pick() to get the current atom under the cursor.
this.parent.addEventListener('click', function(event) {
     	console.log('inside');
	 this.viewer.forEach(function(go) {
    go.colorBy(this.color.uniform('red'), go.select({rindices : [1,2,3,4,5,6,]}));
     });
  this.viewer.requestRedraw();
    var rect = this.viewer.boundingClientRect();
    var picked = this.viewer.pick({ x : event.clientX - rect.left,
                               y : event.clientY - rect.top });
    if (prevPicked !== null && picked !== null &&
        picked.target() === prevPicked.atom) {
      return;
    }
    if (prevPicked !== null) {
      // reset color of previously picked atom.
      setColorForAtom(prevPicked.node, prevPicked.atom, prevPicked.color);
    }
    if (picked !== null) {
      var atom = picked.target();
      var atoms = atom.residue()._r;
      console.log(atoms);
      document.getElementById('picked-atom-name').innerHTML = atom.qualifiedName();
      // get RGBA color and store in the color array, so we know what it was
      // before changing it to the highlight color.
      var color = [0,0,0,0];
      picked.node().getColorForAtom(atom, color);
      prevPicked = { atom : atom, color : color, node : picked.node() };

      setColorForAtom(picked.node(), atoms, 'red');
    } else {
      document.getElementById('picked-atom-name').innerHTML = '&nbsp;';
      prevPicked = null;
    }
    this.viewer.requestRedraw();
});*/


}
