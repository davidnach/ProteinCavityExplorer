


declare var pv:any;


import { Component, OnInit } from '@angular/core';

@Component({
  selector: 'app-experiment-results',
  templateUrl: './experiment-results.component.html',
  styleUrls: ['./experiment-results.component.css']
})
export class ExperimentResultsComponent implements OnInit {




  constructor() {

  }


  ngOnInit() {
        var options = {
          width: 600,
          height: 600,
          antialias: true,
          quality : 'high'
        };
          var viewer = pv.Viewer(document.getElementById("viewer"), options);
          var url = "http://www.rcsb.org/pdb/files/1ake.pdb";
          pv.io.fetchPdb(url, function(structure) {
            viewer.cartoon('protein', structure);
            viewer.autoZoom();
            viewer.centerOn(structure);
          });
 }




}
