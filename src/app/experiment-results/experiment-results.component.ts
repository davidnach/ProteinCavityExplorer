
declare var pv:any;

import {UserInputService} from '../userinput.service';
import { Component, OnInit, } from '@angular/core';
import { GlobalData} from '../globaldata.service';

@Component({
  selector: 'app-experiment-results',
  templateUrl: './experiment-results.component.html',
  styleUrls: ['./experiment-results.component.css']
})
export class ExperimentResultsComponent implements OnInit {
    pid : string;
    experimentId : string;
    
  constructor(private userInputService : UserInputService,
private globaldata : GlobalData) {

  }


  ngOnInit() {
	  this.experimentId = this.globaldata.getExperimentId();
          this.pid = this.globaldata.getPid();
          console.log(this.pid)
          var options = {
          width: 600,
          height: 600,
          antialias: true,
          quality : 'high'
        };
	  var parent = document.getElementById('viewer');
          var viewer = pv.Viewer(parent, options);
          var url = "http://www.rcsb.org/pdb/files/";
          url = url + this.pid + ".pdb";
          pv.io.fetchPdb(url, function(structure) {
            viewer.cartoon('protein', structure);
            viewer.autoZoom();
            viewer.centerOn(structure);

	 });
	 }

}
