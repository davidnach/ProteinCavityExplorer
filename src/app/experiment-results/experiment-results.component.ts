


declare var pv:any;

import {UserInputService} from '../userinput.service';
import { Component, OnInit } from '@angular/core';

@Component({
  selector: 'app-experiment-results',
  templateUrl: './experiment-results.component.html',
  styleUrls: ['./experiment-results.component.css']
})
export class ExperimentResultsComponent implements OnInit {
    pid : string;



  constructor(private userInputService : UserInputService) {

  }


  ngOnInit() {
          this.pid = this.userInputService.getPid();
          console.log(this.pid)
          var options = {
          width: 600,
          height: 600,
          antialias: true,
          quality : 'high'
        };
          var viewer = pv.Viewer(document.getElementById("viewer"), options);
          var url = "http://www.rcsb.org/pdb/files/";
          url = url + this.pid + ".pdb";
          pv.io.fetchPdb(url, function(structure) {
            viewer.cartoon('protein', structure);
            viewer.autoZoom();
            viewer.centerOn(structure);
          });
 }




}
