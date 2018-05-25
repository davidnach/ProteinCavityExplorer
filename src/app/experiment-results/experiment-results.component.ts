
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

@Component({
  selector: 'app-experiment-results',
  templateUrl: './experiment-results.component.html',
  styleUrls: ['./experiment-results.component.css']
})
export class ExperimentResultsComponent   {
    pid : string;
    experimentId : string;
    showSpinner: boolean = true;
    residueTypeCount: Map<string, number>;
    tempResidueTypeCount: Map<string,number>;
    color;
    viewer;
    parent;
    structure;
    orig_structure;
    pocketNums;
    pocketNumber:number;
    pocketSelected : boolean = false;


   //d3 stuff
    radius = 10;

  constructor(private userInputService : UserInputService,
	      private globaldata : GlobalData,
	      private http: HttpClient,
	      @Inject(DOCUMENT) document) {
	this.residueTypeCount = new Map<string,number>();
        this.tempResidueTypeCount = new Map<string,number>();
        this.initializeTempResidueTypeCount();
  }


  ngOnInit() {
	this.pocketNums = Array.from(this.globaldata.getMap().keys());
	this.experimentId = this.globaldata.getExperimentId();
	console.log(this.experimentId);
	this.pid = this.globaldata.getPid();
	console.log(this.pid)
	this.parent = document.getElementById('viewer');
	this.color = pv.color;

    
	var options = {
          width: 600,
          height: 600,
          antialias: true,
          quality : 'high'
         };
 	
	this.viewer = pv.Viewer(this.parent, options);
 	this.color = pv.color;
        var url = "http://www.rcsb.org/pdb/files/";
        url = url + this.pid + ".pdb";
        pv.io.fetchPdb(url, (structure) => {
	    this.orig_structure = structure;
            //this.viewer.cartoon('protein', structure);
            this.viewer.autoZoom();
            this.viewer.centerOn(structure);

	});
      // print PDB into box
      //var fileinfo = document.getElementById("fileinfo");
      //fileinfo.innerHTML(this.http.get(url));
	
}


 colorPocketResidues(pocketNum : number){
     this.viewer.clear();
    this.structure = this.orig_structure;
   // this.viewer.cartoon('protein',this.structure);
    var residues = this.globaldata.getPocketResidues(pocketNum);
    var sel = this.structure.select({rindices: residues});
    this.fillResidueTypeCount(sel);
    this.viewer.spheres('structure.sel', sel,{}); 
    this.viewer.forEach( (structure) => {
	structure.colorBy(this.color.uniform('red'), sel);
    });
    //this.viewer.requestRedraw();
    this.pocketSelected = true;
    
}


fillResidueTypeCount(residueSubset){
    this.tempResidueTypeCount = new Map<string,number>();
    this.initializeTempResidueTypeCount();
    residueSubset.eachResidue( (r) => {
		console.log(r.name());
   		this.tempResidueTypeCount.set(r.name(),this.tempResidueTypeCount.get(r.name()) + 1);
    });
    this.residueTypeCount = this.tempResidueTypeCount;
 
}

initializeTempResidueTypeCount(){
    this.tempResidueTypeCount.set("ALA",0);
    this.tempResidueTypeCount.set("GLY",0);
    this.tempResidueTypeCount.set("ILE",0);
    this.tempResidueTypeCount.set("LEU",0);
    this.tempResidueTypeCount.set("PRO",0);
    this.tempResidueTypeCount.set("VAL",0);
    this.tempResidueTypeCount.set("PHE",0);
    this.tempResidueTypeCount.set("TRP",0);
    this.tempResidueTypeCount.set("TYR",0);
    this.tempResidueTypeCount.set("ASP",0);
    this.tempResidueTypeCount.set("GLU",0);
    this.tempResidueTypeCount.set("ARG",0);
    this.tempResidueTypeCount.set("HIS",0);
    this.tempResidueTypeCount.set("LYS",0);
    this.tempResidueTypeCount.set("SER",0);
    this.tempResidueTypeCount.set("THR",0);
    this.tempResidueTypeCount.set("CYS",0);
    this.tempResidueTypeCount.set("MET",0);
    this.tempResidueTypeCount.set("ASN",0);
    this.tempResidueTypeCount.set("GLN",0);
    
}

}

