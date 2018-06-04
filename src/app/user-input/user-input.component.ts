import { Component, OnInit, ElementRef,EventEmitter,Output} from '@angular/core';
import {NgForm,AbstractControl,ValidatorFn} from '@angular/forms';
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
  pidSet: Set<string>;
  pidExists: boolean = false;
  submittedForm: boolean = false;
  constructor(private userInputService : UserInputService,
	      private http: HttpClient,
	      private serverService : ServerService,
	      private globalData : GlobalData) {
      console.log('user-input.component.ts constructor initialized.');
      
      var url = 'https://www.rcsb.org';
      var opts = '/pdb/rest/getCurrent';
      let parseString = require('xml2js').parseString;
      var xml;
      this.pidSet = new Set<string>();
      this.http.get(url+opts, {responseType:'text'})
	  .map(response => {
	      let data;
	      parseString(response, (err, result) => {
		  console.dir(result); // print JSON object
		  var i;
	          for(i = 0; i < result.current['PDB'].length; i++){
			this.pidSet.add(result.current['PDB'][i]['$']['structureId']);
		  }
		  
		  data = result;
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
	this.checkIfPidExists(form.value.pid);
	this.submittedForm = true;
	//form.value.pid = form.value.pid.toLowerCase();
	console.log('user entered ' + form.value.pid + 'from user component');
	if (form.valid && this.pidExists ){ 
	    this.userInputService.pidEntered.emit(form.value.pid.toLowerCase());
	    this.globalData.setPid(form.value.pid);
	}
	form.resetForm();
    }
    
    submitExperimentId(form : NgForm){
         this.submittedForm=true;
	 if(form.valid){
		this.userInputService.experimentIdEntered.emit(form.value.experimentId);
         	this.globalData.setExperimentId(form.value.experimentId);
    	 }
	 
}

    checkIfPidExists(pid){
	console.log(pid.toUpperCase());
	this.pidExists = this.pidSet.has(pid.toUpperCase());
	console.log(this.pidExists);	
    }	

}
