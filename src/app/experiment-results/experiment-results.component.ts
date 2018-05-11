
declare var pv:any;
declare var structure:any;
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
    showSpinner: boolean = true;
    color;
    viewer;
    parent;
    pocketNums;
    pocketNumber:number;


  constructor(private userInputService : UserInputService,
	      private globaldata : GlobalData) {
 	this.pocketNums = Array.from(this.globaldata.getMap().keys());
	this.experimentId = this.globaldata.getExperimentId();
	console.log(this.pocketNums);
  }


  ngOnInit() {
    console.log(this.experimentId);
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

 colorPocketResidues(pocketNum : number){
    var residues = this.globaldata.getPocketResidues(pocketNum);
    var color_ = this.color;
    this.viewer.forEach(function(go) {
	console.log(go);
	go.colorBy(color_.uniform('red'), go.select({rindices : residues}));
    });
    this.viewer.requestRedraw();
}

}
