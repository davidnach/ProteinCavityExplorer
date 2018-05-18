
declare var pv:any;
declare var structure:any;
import {UserInputService} from '../userinput.service';
import { Component, OnInit, } from '@angular/core';
import { GlobalData} from '../globaldata.service';
import { HttpClient, HttpErrorResponse} from '@angular/common/http';
import { Inject }  from '@angular/core';
import { DOCUMENT } from '@angular/common';

@Component({
  selector: 'app-experiment-results',
  templateUrl: './experiment-results.component.html',
  styleUrls: ['./experiment-results.component.css']
})
export class ExperimentResultsComponent   {
    pid : string;
    experimentId : string;
    showSpinner: boolean = true;
    color;
    viewer;
    parent;
    structure;
    orig_structure;
    pocketNums;
    pocketNumber:number;


  constructor(private userInputService : UserInputService,
	      private globaldata : GlobalData,
	      private http: HttpClient,
	      @Inject(DOCUMENT) document) {
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
            this.viewer.cartoon('protein', structure);
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
    this.viewer.cartoon('protein',this.structure);
    var residues = this.globaldata.getPocketResidues(pocketNum);
    var sel = this.structure.select({rindices: residues});
    this.viewer.spheres('structure.sel', sel,{}); 
    this.viewer.forEach( (structure) => {
	structure.colorBy(this.color.uniform('red'), sel);
    });
    this.viewer.requestRedraw();
}

}


