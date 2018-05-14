import { Component, OnInit, ElementRef,EventEmitter,Output} from '@angular/core';
import {NgForm} from '@angular/forms';
import { HttpClient , HttpErrorResponse} from '@angular/common/http';
import {UserInputService} from '../userinput.service';
import {ServerService} from '../server.service';
import {GlobalData} from '../globaldata.service';

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
	}// else {
	// validation failed
	//}
	//need make sure pid exists
	// if(form.valid){
	//   this.serverService.checkIfPidExists("1ake").subscribe(
	//     (error : HttpErrorResponse) => {console.log(error.status);
	//     console.log("something");})
	//   }
    }
/*
    checkPidValidity(pid : string) {
	console.log('checking pid validity...');
	var url = 'https://www.rcsb.org/pdb/rest/idStatus?structureId=' + pid + '';
	var body = '';
	this.http.post(url, body, {
	    headers: new HttpHeaders()
	        .set('Content-Type', 'text/xml')
	        .append('Access-Control-Allow-Methods', 'GET,PUT,POST,DELETE,PATCH,OPTIONS')
	        .append('Access-Control-Allow-Origin', '*')
	        .append('Access-Control-Allow-Headers', "Access-Control-Allow-Headers, Access-Control-Allow-Origin, Access-Control-Request-Method")
	    , responseType:'text'
	}).subscribe(data => {
	    console.log(data);
	},
	(err : HttpErrorResponse) => {
	    if (err.error instanceOf Error) {
		// client side or network error
		console.log('Error: ', err.error.message);
	    } else {
		// backend returned unsuccessful response code
		console.log(`Backend returned code ${err.status}, body was: ${err.error}`);
	    }
	}
      );
    }*/
      //this.serverService.getPidStatus(form.value.pid).subscribe(
      // data => { this.foods = data },
    submitExperimentId(form : NgForm){
	if(form.valid) this.userInputService.experimentIdEntered.emit(form.value.experimentId);
	this.globalData.setExperimentId(form.value.experimentId);
    }
}
