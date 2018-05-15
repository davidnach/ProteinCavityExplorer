import { Component, OnInit, ElementRef,EventEmitter,Output} from '@angular/core';
import {NgForm} from '@angular/forms';
import { HttpClient , HttpErrorResponse} from '@angular/common/http';
import {UserInputService} from '../userinput.service';
import {ServerService} from '../server.service';
import {GlobalData} from '../globaldata.service';
import { Observable } from 'rxjs/Rx';
import 'rxjs/add/operator/map';

@Component({
  selector: 'app-user-input',
  templateUrl: './user-input.component.html',
  styleUrls: ['./user-input.component.css']
})
export class UserInputComponent implements OnInit {

  constructor(private userInputService : UserInputService,
	      private http: HttpClient,
	      private serverService : ServerService,
	      private globalData : GlobalData) {
      console.log('user-input.component.ts constructor initialized.');
      
      var url = 'https://www.rcsb.org';
      var opts = '/pdb/rest/getCurrent';
      let parseString = require('xml2js').parseString;
      var xml;
      this.http.get(url+opts, {responseType:'text'})
	  .map(response => {
	      let data;
	      parseString(response, function(err, result) {
		  console.dir(result); // print JSON object
		  data = result
	      });
	      return data;
	  })
	  .subscribe(data => {
	      console.log(data); // print again...
	  }); 
  }
    
    ngOnInit() {
	console.log('user-input.component onInit opened');
    }

    submitPid(form : NgForm){
	//form.value.pid = form.value.pid.toLowerCase();
	console.log('user entered ' + form.value.pid + 'from user component');
	if (form.valid){// && checkPidValidity(form.value.pid)) {
	    this.userInputService.pidEntered.emit(form.value.pid);
	    this.globalData.setPid(form.value.pid);
	}
    }
    
    submitExperimentId(form : NgForm){
	if(form.valid) this.userInputService.experimentIdEntered.emit(form.value.experimentId);
	this.globalData.setExperimentId(form.value.experimentId);
    }

}
