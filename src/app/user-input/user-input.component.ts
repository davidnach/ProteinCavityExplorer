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
  private http: HttpClient, private serverService : ServerService,
private globalData : GlobalData) { }

  ngOnInit() {
  }

  submitPid(form : NgForm){
    console.log('user entered ' + form.value.pid + 'from user component');
    if(form.valid){
      this.userInputService.pidEntered.emit(form.value.pid);
      this.globalData.setPid(form.value.pid);
    }
      //need make sure pid exists
      // if(form.valid){
      //   this.serverService.checkIfPidExists("1ake").subscribe(
      //     (error : HttpErrorResponse) => {console.log(error.status);
      //     console.log("something");})
      //   }
    }


  submitExperimentId(form : NgForm){
      if(form.valid) this.userInputService.experimentIdEntered.emit(form.value.experimentId);
      this.globalData.setExperimentId(form.value.experimentId);
  }
}
