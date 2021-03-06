
declare var pv:any;
declare var structure:any;

import * as d3 from 'd3';
import {UserInputService} from '../userinput.service';
import { Component, OnInit, } from '@angular/core';
import { GlobalData} from '../globaldata.service';
import { HttpClient, HttpErrorResponse} from '@angular/common/http';
import { Inject }  from '@angular/core';
import { DOCUMENT } from '@angular/common';
import { BarchartsComponent} from '../barcharts/barcharts.component';
import { PieChartComponent} from '../piechart/piechart.component';
import {MatCheckboxModule} from '@angular/material/checkbox';
import {VERSION} from '@angular/material';

@Component({
  selector: 'app-experiment-results',
  templateUrl: './experiment-results.component.html',
  styleUrls: ['./experiment-results.component.css']  
})
export class ExperimentResultsComponent   {
    aminos = ['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL'];
    version = VERSION;
    isLinear = true;
  items = [
  ];
    pid : string;
    experimentId : string;
    showSpinner: boolean = true;
    residueTypeCount: Map<string, number>;
    baseResidueTypeCount: Map<string,number>;
    color;
    viewer;
    parent;
    structure;
    orig_structure;
    selectedPocketsTracker = [];
    pocketNums = [];
    pocketSelected : boolean = false;
    barChartData : any[] = [];
    selectedResidues;
    numPockets = 0;
   //d3 stuff
    radius = 10;
    initGraphSelect = 0;
    selectedPocketArea;
    pdburl;

  constructor(private userInputService : UserInputService,
	      private globaldata : GlobalData,
	      private http: HttpClient,
	      @Inject(DOCUMENT) document) {
	      this.residueTypeCount = new Map <string,number>();
  }


  ngOnInit() {
     	  
        this.setPocketNums();
        this.experimentId = this.globaldata.getExperimentId();
     	
       console.log(this.experimentId);
	this.pid = this.globaldata.getPid();
	this.pdburl = "https://files.rcsb.org/view/" + this.pid + ".pdb";
	console.log(this.pdburl);
	console.log(this.pid)
	this.parent = document.getElementById('viewer');
	this.color = pv.color;
    
	var options = {
          width: 500,
          height: 500,
          antialias: true,
          quality : 'high'
         };
 	
	this.viewer = pv.Viewer(this.parent, options);
 	this.color = pv.color;
        
        var url = "http://www.rcsb.org/pdb/files/";
        url = url + this.pid + ".pdb";
        pv.io.fetchPdb(url, (structure) => {
	    this.orig_structure = structure;
            this.viewer.cartoon('protein', structure);
            this.viewer.autoZoom();
            this.viewer.centerOn(structure);
	    this.changed();
	    //this.selectAll(); // display fully selected graph on page load
	});

var prevPicked = null;      

this.parent.addEventListener('mousemove', (event) => {
    var rect = this.viewer.boundingClientRect();
    var picked = this.viewer.pick({ x : event.clientX - rect.left,
                               y : event.clientY - rect.top });
    if (prevPicked !== null && picked !== null &&
        picked.target() === prevPicked.atom) {
      return;
    }
    if (prevPicked !== null) {
      // reset color of previously picked atom.
      this.setColorForAtom(prevPicked.node, prevPicked.atom, prevPicked.color);
    }
    if (picked !== null) {
      var atom = picked.target();
      document.getElementById('picked-atom-name').innerHTML = atom.qualifiedName();
      // get RGBA color and store in the color array, so we know what it was
      // before changing it to the highlight color.
      var color = [0,0,0,0];
      picked.node().getColorForAtom(atom, color);
      prevPicked = { atom : atom, color : color, node : picked.node() };

      this.setColorForAtom(picked.node(), atom, 'blue');
    } else {
      document.getElementById('picked-atom-name').innerHTML = '&nbsp;';
      prevPicked = null;
    }
   this.viewer.requestRedraw();
});
  }
  
   

    setColorForAtom(go, atom, color) {
    var view = go.structure().createEmptyView();
    view.addAtom(atom);
    go.colorBy(pv.color.uniform(color), view);
   }
    
    colorPocketResidues(pocketNum : number){
     this.selectedPocketArea = this.globaldata.getPocketArea(pocketNum);
     this.viewer.clear();
     this.structure = this.orig_structure;
     this.viewer.cartoon('protein',this.structure);
     var residues = this.globaldata.getPocketResidues(pocketNum);
     var sel = this.structure.select({rindices: residues});
     //this.fillResidueTypeCount(sel);
     this.viewer.spheres('structure.sel', sel,{});
     this.viewer.forEach( (structure) => {
	 structure.colorBy(this.color.uniform('red'), sel);
     });
    this.structure = this.orig_structure;
    this.viewer.cartoon('protein',this.structure);
   /* this.selectedResidues = [];
    var i;
    var j;
    for(i = 0; i < this.pocketNums; i++){
	if(this.items[i].checked == true){
		var pocketResidues = this.globaldata.getPocketResidues(i);
		for(j = 0; j < pocketResidues.length;j++){
			this.selectedResidues.push(pocketResidues[i]);
		}
        }
    }
    var sel = this.structure.select({rindices: this.selectedResidues});
    this.fillResidueTypeCount(sel);*/
    this.viewer.spheres('structure.sel', sel,{}); 
    this.viewer.forEach( (structure) => {
	structure.colorBy(this.color.uniform('red'), sel);
    }); 
    this.pocketSelected = true;
 }


setPocketNums(){
  this.pocketNums = Array.from(this.globaldata.getMap().keys());
  this.numPockets = this.pocketNums.length;
  var i;
  for(i = 0; i < this.numPockets; i++){
	this.selectedPocketsTracker.push({checked : false});
  } 
  this.selectedPocketsTracker[0].checked = true;
}


fillResidueTypeCount(residueSubset){
    this.residueTypeCount = new Map<string, number>();
    this.initializeBaseResidueTypeCount(this.residueTypeCount);
    var residueCount = 0;
    console.log(residueSubset);
    console.log(residueSubset.residueCount());
    residueSubset.eachResidue( (r) => {
		residueCount++;
		console.log(r.name());
   		this.residueTypeCount.set(r.name(),this.residueTypeCount.get(r.name()) + 1);
    });
    console.log(residueCount);
}

initializeBaseResidueTypeCount(resTypeCount){
    resTypeCount.set("ALA",0);
    resTypeCount.set("GLY",0);
    resTypeCount.set("ILE",0);
    resTypeCount.set("LEU",0);
    resTypeCount.set("PRO",0);
    resTypeCount.set("VAL",0);
    resTypeCount.set("PHE",0);
    resTypeCount.set("TRP",0);
    resTypeCount.set("TYR",0);
    resTypeCount.set("ASP",0);
    resTypeCount.set("GLU",0);
    resTypeCount.set("ARG",0);
    resTypeCount.set("HIS",0);
    resTypeCount.set("LYS",0);
    resTypeCount.set("SER",0);
    resTypeCount.set("THR",0);
    resTypeCount.set("CYS",0);
    resTypeCount.set("MET",0);
    resTypeCount.set("ASN",0);
    resTypeCount.set("GLN",0);
}

selectAll(){
    console.log('selecting all...struct: ' + this.orig_structure);
    var i;
    for(i = 0; i < this.numPockets; i++) this.selectedPocketsTracker[i].checked = true;
    this.changed();
}

unselectAll(){
   var i;
   for(i = 0; i < this.numPockets; i++) this.selectedPocketsTracker[i].checked = false;
   this.changed();
}

changed(){
    this.structure = this.orig_structure;
    console.log(this.structure);
    var pocketResidueTypeCounts = [];
    var loadBarChartData = [];
    var i;
    var j;
    var loadDataTracker = 0;
    var residues = [];
    for(i = 0; i < this.numPockets; i++){
	if(this.selectedPocketsTracker[i].checked == true){
		residues = this.globaldata.getPocketResidues(i);
 		console.log(residues);
		var sel = this.structure.select({rindices: residues});
		this.fillResidueTypeCount(sel);
                console.log(this.residueTypeCount);
		pocketResidueTypeCounts = [];
                for(j = 0; j < this.aminos.length; j++){
			pocketResidueTypeCounts[j] = this.residueTypeCount.get(this.aminos[j]);
		}
		loadBarChartData[loadDataTracker]= {data : pocketResidueTypeCounts, label: 'pocket' + i};		
		loadDataTracker++;
                residues = [];	
	}
    }
    this.barChartData = loadBarChartData;
    console.log(this.barChartData);		 
}

}

